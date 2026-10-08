 module smbpal

    use nml
    use smbpal_precision
    use insolation
    use interp_time 
    use ncio
    use smb_pdd 
    use smb_itm 

    implicit none 

    integer, parameter :: ndays = 360   ! 360-day year
    integer, parameter :: ndays_mon = 30   ! 30 days per month  

    ! itm: daily forcing interpolated from the monthly fields at ndays_daily
    ! days, snowpack budget every dt_itm days, annual PDDs from every dt_pdd days
    integer, parameter :: ndays_daily = 37
    integer, parameter :: dt_itm      = 2
    integer, parameter :: nstep_itm   = ndays/dt_itm
    integer, parameter :: dt_pdd      = 10
    integer, parameter :: nstep_pdd   = ndays/dt_pdd
        
    type smbpal_param_class
        type(itm_par_class) :: itm
        logical    :: const_insol
        real(prec) :: const_kabp
        character(len=512)  :: insol_fldr 
        character(len=16)   :: abl_method 
        real(prec) :: sigma_snow, sigma_melt, sigma_land
        real(prec) :: sf_a, sf_b, firn_fac  
        real(prec) :: mm_snow, mm_ice 

        real(prec), allocatable :: x(:), y(:)
        real(prec), allocatable :: lats(:,:)           ! Latitude of domain [deg N]

    end type 

    type smbpal_state_class 
        real(prec), allocatable   :: t2m(:,:)            ! Surface temperature [K]
        real(prec), allocatable   :: pr(:,:), sf(:,:)    ! Precip, snowfall [mm/a or mm/d]
        real(prec), allocatable   :: S(:,:)              ! Insolation [W/m2]
        real(prec), allocatable   :: sigma(:,:)          ! Effective temp. (ie, PDDs) [num. of days]
        real(prec), allocatable   :: PDDs(:,:)           ! Effective temp. (ie, PDDs) [num. of days]
        real(prec), allocatable   :: tsrf(:,:)           ! Effective temp. (ie, PDDs) [num. of days]
        
        ! Prognostic variables
        real(prec), allocatable   :: H_snow(:,:)         ! Snow thickness [mm]
        real(prec), allocatable   :: alb_s(:,:)          ! Surface albedo 
        real(prec), allocatable   :: smbi(:,:), smb(:,:) ! Surface mass balance [mm/a or mm/d]
        real(prec), allocatable   :: melt(:,:), runoff(:,:), refrz(:,:)   ! smb components
        real(prec), allocatable   :: melt_net(:,:)       ! Net surface melt, for calculating surface temp [mm]
    end type 

    type smbpal_class
        type(smbpal_param_class) :: par 
        type(smbpal_state_class) :: now, mon(12), ann
    end type

    type smbpal_point_class
        ! smbpal state at one grid point (the fields of smbpal_state_class)
        real(prec) :: t2m = 0.0, pr = 0.0, sf = 0.0, S = 0.0
        real(prec) :: sigma = 0.0, PDDs = 0.0, tsrf = 0.0
        real(prec) :: H_snow = 0.0, alb_s = 0.0, smbi = 0.0, smb = 0.0
        real(prec) :: melt = 0.0, runoff = 0.0, refrz = 0.0, melt_net = 0.0
    end type 

    type smbpal_itm_forcing_class
        ! Daily forcing of one itm year (smbpal_itm_forcing_init)
        integer :: days(ndays_daily)                ! Days of the daily forcing fields
        integer :: k_step(nstep_itm)                ! Forcing day index following each itm step
        integer :: k_pdd(nstep_pdd)                 ! Forcing day index following each pdd step
        real(prec), allocatable :: t2m(:,:,:)       ! [K] Daily temperature [nx,ny,ndays_daily]
        real(prec), allocatable :: pr(:,:,:)        ! [mm we/d] Daily precipitation
        real(prec), allocatable :: sf(:,:,:)        ! [mm we/d] Daily snowfall
        type(insol_day_class)   :: insol(nstep_itm) ! Insolation of each itm step
    end type 

    interface operator(+)
        module procedure smbpal_point_add
    end interface 

    interface operator(/)
        module procedure smbpal_point_div
    end interface 
    
    private
    public :: smbpal_class
    public :: smbpal_init 
    public :: smbpal_update_2temp, smbpal_update_monthly 
    public :: smbpal_update_monthly_equil
    public :: smbpal_end
    public :: smbpal_write_init, smbpal_write
    public :: smbpal_restart_write, smbpal_restart_read

contains 

    subroutine smbpal_init(smb,filename,x,y,lats,group,itm_group)

        implicit none 

        type(smbpal_class) :: smb
        character(len=*), intent(IN)  :: filename  ! Parameter file 
        real(prec) :: x(:), y(:), lats(:,:)
        character(len=*),  intent(IN), optional :: group, itm_group

        ! Local variables
        integer :: nx, ny, m  
        real(prec) :: tmp 
        character(len=32) :: nml_group, itm_nml_group

        ! Make sure we know the namelist group for the smbpal block
        if (present(group)) then
            nml_group = trim(group)
        else
            nml_group = "smbpal"         ! Default parameter blcok name
        end if

        ! Make sure we know the namelist group for the itm block
        if (present(itm_group)) then
            itm_nml_group = trim(itm_group)
        else
            itm_nml_group = "itm"         ! Default parameter blcok name
        end if

        nx = size(x,1)
        ny = size(y,1)

        ! Load smbpal parameters
        call smbpal_par_load(smb%par,filename,init=.TRUE.,group=nml_group,itm_group=itm_nml_group)

        ! Additionally define dimension info 
        if (allocated(smb%par%x)) deallocate(smb%par%x)
        if (allocated(smb%par%y)) deallocate(smb%par%y)
        if (allocated(smb%par%lats)) deallocate(smb%par%lats)
        allocate(smb%par%x(nx),smb%par%y(ny),smb%par%lats(nx,ny))

        smb%par%x    = x 
        smb%par%y    = y 
        smb%par%lats = lats 

        ! Allocate the smbpal object 
        call smbpal_allocate(smb%now,nx,ny)
        call smbpal_allocate(smb%ann,nx,ny)
    
        do m = 1, 12 
            call smbpal_allocate(smb%mon(m),nx,ny)
        end do 

        ! Initialize the state variables 
        smb%now%H_snow = smb%par%itm%H_snow_max 

        ! Test calculation of insolation to load orbital params 
        tmp = calc_insol_day(180,65.d0,0.d0,fldr=smb%par%insol_fldr)

        return 

    end subroutine smbpal_init

    subroutine smbpal_update_2temp(smb,t2m_ann,t2m_sum,pr_ann,z_srf,H_ice,time_bp,sf_ann, &
                                   file_out,file_out_mon,write_init,calc_mon,write_now)
        ! Generate climate using two points in year (Tsum,Tann)

        implicit none 
        
        type(smbpal_class), intent(INOUT) :: smb
        real(prec), intent(IN) :: t2m_ann(:,:), t2m_sum(:,:)
        real(prec), intent(IN) ::  pr_ann(:,:), z_srf(:,:), H_ice(:,:)
        real(prec), intent(IN) :: time_bp       ! years BP 
        real(prec), intent(IN), optional :: sf_ann(:,:)
        character(len=*), intent(IN), optional :: file_out      ! Annual output
        character(len=*), intent(IN), optional :: file_out_mon  ! Monthly output
        logical, intent(IN), optional :: write_init, calc_mon, write_now

        ! Local variables
        real(prec), allocatable :: t2m(:,:,:), pr(:,:,:), sf(:,:,:) 
        integer :: day, m
        real(prec) :: dt 

        allocate(t2m(size(t2m_ann,1),size(t2m_ann,2),12))
        allocate( pr(size(t2m_ann,1),size(t2m_ann,2),12))
        allocate( sf(size(t2m_ann,1),size(t2m_ann,2),12))
        
        do m = 1, 12
            ! Determine t2m, pr, sf and S today 
            day = m*30
            t2m(:,:,m) = t2m_ann-(t2m_sum-t2m_ann)*cos(2.0*pi*real(day-15)/real(ndays))
        
            pr(:,:,m)  = pr_ann
            if (present(sf_ann)) then 
                sf(:,:,m) = sf_ann
            else 
                sf(:,:,m) = pr(:,:,m) * calc_snowfrac(t2m(:,:,m),smb%par%sf_a,smb%par%sf_b)
            end if 

        end do  

        ! Call monthly interface
        call smbpal_update_monthly(smb,t2m,pr,z_srf,H_ice,time_bp,sf, &
                        file_out,file_out_mon,write_init,calc_mon,write_now)
        return 

    end subroutine smbpal_update_2temp

    subroutine smbpal_update_monthly_equil(smb,t2m,pr,z_srf,H_ice,time_bp,time_equil,sf)
        ! Generate climate using monthly input data [nx,ny,nmon]
        
        implicit none 
        
        type(smbpal_class), intent(INOUT) :: smb
        real(prec), intent(IN) :: t2m(:,:,:), pr(:,:,:)
        real(prec), intent(IN) ::  z_srf(:,:), H_ice(:,:)
        real(prec), intent(IN) :: time_bp       ! years BP 
        real(prec), intent(IN) :: time_equil    ! years to equilibrate
        real(prec), intent(IN), optional :: sf(:,:,:)

        ! Local variables
        integer :: n 
        type(smbpal_itm_forcing_class) :: frc 

        ! Loop over equilibration years to update snowpack thickness 
        if (trim(smb%par%abl_method) .eq. "itm") then 
            ! The daily forcing is the same every year: prepare it once
            call smbpal_itm_forcing_init(frc,smb%par,t2m,pr,time_bp,sf)
            do n = 1, int(time_equil) 
                call smbpal_update_itm(smb,frc,z_srf,H_ice,calc_monthly=.FALSE.)
            end do 
        else 
            do n = 1, int(time_equil) 
                call smbpal_update_monthly(smb,t2m,pr,z_srf,H_ice,time_bp,sf)
            end do 
        end if 

        return 

    end subroutine smbpal_update_monthly_equil


    subroutine smbpal_update_monthly(smb,t2m,pr,z_srf,H_ice,time_bp,sf, &
                        file_out,file_out_mon,write_init,calc_mon,write_now)
        ! Generate climate using monthly input data [nx,ny,nmon]
        
        implicit none 
        
        type(smbpal_class), intent(INOUT) :: smb
        real(prec),         intent(IN) :: t2m(:,:,:)                ! [K] Monthly temperature fields
        real(prec),         intent(IN) :: pr(:,:,:)                 ! [mm we/d] Monthly precipitation rate fields 
        real(prec),         intent(IN) :: z_srf(:,:)                ! [m] Surface elevation 
        real(prec),         intent(IN) :: H_ice(:,:)                ! [m] Ice thickness 
        real(prec),         intent(IN) :: time_bp                   ! [years BP] Current time (for insolation) 
        real(prec),         intent(IN), optional :: sf(:,:,:)       ! [mm we/d] Monthly snowfall rate fields 
        character(len=*),   intent(IN), optional :: file_out        ! Annual output filename
        character(len=*),   intent(IN), optional :: file_out_mon    ! Monthly output filename
        logical,            intent(IN), optional :: write_init      ! Flag for whether to initialize writing of output file
        logical,            intent(IN), optional :: calc_mon        ! Flag for whether to calculate monthly averages (itm only)
        logical,            intent(IN), optional :: write_now       ! Flag for whether to write the current time to file

        ! Local variables
        logical :: init_now, write_out_now
        logical :: calc_monthly, init_mon, write_mon
        integer :: k, m 
        type(smbpal_itm_forcing_class) :: frc 
        
        real(prec), allocatable :: tmp4(:,:)
        real(prec), allocatable :: t2m_ann(:,:), pr_ann(:,:), sf_ann(:,:) 
        real(prec), allocatable :: PDDs_ann(:,:) 

        write_out_now = .FALSE. 
        if (present(write_now) .and. present(file_out)) write_out_now = write_now 

        ! Determine whether this is first time running (for output)
        init_now = .FALSE. 
        if (write_out_now .and. present(write_init)) init_now = write_init 

        allocate(tmp4(size(t2m,1),size(t2m,2)))

        if (trim(smb%par%abl_method) .eq. "itm") then 
            
            calc_monthly = .FALSE. 
            if (present(calc_mon))     calc_monthly = calc_mon 
            if (present(file_out_mon)) calc_monthly = .TRUE. 

            ! Generate daily climate from monthly input and run the year
            call smbpal_itm_forcing_init(frc,smb%par,t2m,pr,time_bp,sf)
            call smbpal_update_itm(smb,frc,z_srf,H_ice,calc_monthly)

            ! Monthly I/O 
            write_mon = .FALSE. 
            if (present(write_now)) write_mon = write_now 
            init_mon = .FALSE. 
            if (present(write_init)) init_mon = write_init 

            if (write_mon .and. calc_monthly .and. present(file_out_mon)) then
                if (init_mon) call smbpal_write_init(smb%par,file_out_mon,z_srf,H_ice)
                do m = 1, 12
                    call smbpal_write(smb%mon(m),file_out_mon,time_bp=time_bp,step="mon",nstep=m)
                end do 
            end if 

        else
            ! PDD method 

            allocate(t2m_ann(size(t2m,1),size(t2m,2)))
            allocate(pr_ann(size(t2m,1),size(t2m,2)))
            allocate(sf_ann(size(t2m,1),size(t2m,2)))
            allocate(PDDs_ann(size(t2m,1),size(t2m,2)))

            t2m_ann = sum(t2m,dim=3) / 12.0 
            pr_ann  = sum(pr, dim=3) / 12.0 *real(ndays,prec)       ! [mm we/d] => [mm we/a]

            if (present(sf)) then 
                sf_ann  = sum(sf,dim=3) / 12.0 *real(ndays,prec)    ! [mm we/d] => [mm we/a]
            else 
                sf_ann  = pr_ann    ! Should be improved in the future 
            end if 

            ! First calculate PDDs for the whole year (input to pdd)
            PDDs_ann = 0.0 
            do k = 1, 12
                smb%now%t2m  = t2m(:,:,k)

                smb%now%sigma = smb%par%sigma_snow 
                where (z_srf .gt. 0.0 .and. H_ice .eq. 0.0)          smb%now%sigma = smb%par%sigma_land 
                where (H_ice .gt. 0.0 .and. smb%now%t2m .ge. 273.15) smb%now%sigma = smb%par%sigma_melt
                
                call calc_temp_effective(tmp4,smb%now%t2m-273.15,smb%now%sigma)
                PDDs_ann = PDDs_ann + tmp4*30.0

            end do 

            ! Populate the ann object with the now object, then calculate the annual values 
            smb%ann = smb%now 
             
            call smbpal_update_pdd(smb%ann,smb%par,PDDs_ann,z_srf,H_ice,t2m_ann,pr_ann,sf_ann)

            ! Note: annual values are output with units of [mm/a]

        end if 
        
        ! Annual I/O 
        if (write_out_now) then
            if (init_now) call smbpal_write_init(smb%par,file_out,z_srf,H_ice)
            call smbpal_write(smb%ann,file_out,time_bp=time_bp,step="ann")

        end if 

        return 

    end subroutine smbpal_update_monthly

    subroutine smbpal_itm_forcing_init(frc,par,t2m,pr,time_bp,sf)
        ! Daily forcing of one itm year: the monthly fields interpolated to
        ! the forcing days and the insolation of each itm step. It depends
        ! only on the climate, so it is shared by all grid points and years.

        implicit none 

        type(smbpal_itm_forcing_class), intent(OUT) :: frc 
        type(smbpal_param_class),       intent(IN)  :: par 
        real(prec), intent(IN) :: t2m(:,:,:)                ! [K] Monthly temperature fields
        real(prec), intent(IN) :: pr(:,:,:)                 ! [mm we/d] Monthly precipitation rate fields 
        real(prec), intent(IN) :: time_bp                   ! [years BP] Current time (for insolation) 
        real(prec), intent(IN), optional :: sf(:,:,:)       ! [mm we/d] Monthly snowfall rate fields 

        ! Local variables
        integer :: nx, ny, k, n, day 
        real(8) :: insol_time
        double precision, allocatable :: tmp(:,:,:)

        nx = size(t2m,1)
        ny = size(t2m,2)

        ! Define daily days 
        do k = 1, ndays_daily-1 
            frc%days(k) = 1 + (k-1)*(ndays / (ndays_daily-1))
        end do 
        frc%days(ndays_daily) = ndays 

        ! Index of the forcing day following each itm and pdd step
        do n = 1, nstep_itm 
            frc%k_step(n) = idx_today(frc%days,1 + (n-1)*dt_itm)
        end do 
        do n = 1, nstep_pdd 
            frc%k_pdd(n) = idx_today(frc%days,1 + (n-1)*dt_pdd)
        end do 

        ! Generate daily climate from monthly input 
        allocate(tmp(nx,ny,ndays_daily))
        allocate(frc%t2m(nx,ny,ndays_daily))
        allocate(frc%pr(nx,ny,ndays_daily))
        allocate(frc%sf(nx,ny,ndays_daily))

        call convert_monthly_daily_3D(dble(t2m),tmp,days=frc%days)
        frc%t2m = tmp 
        call convert_monthly_daily_3D(dble(pr),tmp,days=frc%days)
        frc%pr = tmp 

        if (present(sf)) then 
            call convert_monthly_daily_3D(dble(sf),tmp,days=frc%days)
            frc%sf = tmp 
            where(frc%sf .lt. 0.0) frc%sf = 0.0 
        else 
            do k = 1, ndays_daily 
                frc%sf(:,:,k) = frc%pr(:,:,k) * calc_snowfrac(frc%t2m(:,:,k),par%sf_a,par%sf_b)
            end do 
        end if 

        ! Determine year to use for insolation calcs
        insol_time = time_bp
        if (par%const_insol) insol_time = par%const_kabp*1e3
        
        ! Insolation of each itm step at the latitude nodes 
        do n = 1, nstep_itm 
            day = 1 + (n-1)*dt_itm
            frc%insol(n) = calc_insol_day_spline(day,insol_time,fldr=par%insol_fldr)
        end do 

        return 

    end subroutine smbpal_itm_forcing_init

    subroutine smbpal_update_itm(smb,frc,z_srf,H_ice,calc_monthly)
        ! One itm year from the daily forcing frc. Grid points are
        ! independent, so the year is integrated point by point in parallel.

        implicit none 
        
        type(smbpal_class),             intent(INOUT) :: smb
        type(smbpal_itm_forcing_class), intent(IN)    :: frc 
        real(prec),                     intent(IN)    :: z_srf(:,:), H_ice(:,:)
        logical,                        intent(IN)    :: calc_monthly

        ! Local variables
        integer :: i, j, m, nx, ny 
        type(smbpal_point_class) :: now, ann, mon(12)

        nx = size(z_srf,1)
        ny = size(z_srf,2)

        !$omp parallel do collapse(2) private(i,j,m,now,ann,mon)
        do j = 1, ny 
        do i = 1, nx 
            
            now = smbpal_point_get(smb%now,i,j)

            call smbpal_itm_point(now,ann,mon,smb%par,frc,i,j,z_srf(i,j),H_ice(i,j),calc_monthly)

            call smbpal_point_set(smb%now,i,j,now)
            call smbpal_point_set(smb%ann,i,j,ann)
            if (calc_monthly) then 
                do m = 1, 12 
                    call smbpal_point_set(smb%mon(m),i,j,mon(m))
                end do 
            end if 

        end do 
        end do 
        !$omp end parallel do

        return 

    end subroutine smbpal_update_itm

    subroutine smbpal_itm_point(now,ann,mon,par,frc,i,j,z_srf,H_ice,calc_monthly)
        ! One itm year at grid point (i,j): the annual PDDs, then the
        ! snowpack budget every dt_itm days, averaged over the year
        ! (and the months, if calc_monthly). 

        implicit none 

        type(smbpal_point_class),       intent(INOUT) :: now        ! State, carried between years
        type(smbpal_point_class),       intent(INOUT) :: ann        ! Annual mean
        type(smbpal_point_class),       intent(INOUT) :: mon(12)    ! Monthly means (if calc_monthly)
        type(smbpal_param_class),       intent(IN)    :: par 
        type(smbpal_itm_forcing_class), intent(IN)    :: frc 
        integer,                        intent(IN)    :: i, j 
        real(prec),                     intent(IN)    :: z_srf, H_ice 
        logical,                        intent(IN)    :: calc_monthly

        ! Local variables
        integer :: n, day, k1, mnow, mday 
        real(prec) :: dt    ! [days]
        real(prec) :: teff 
        
        dt = real(dt_itm,prec)

        ! Set sigma to snow sigma for pdd calcs
        now%sigma = par%sigma_snow

        ! First calculate PDDs for the whole year (input to itm)
        now%PDDs = 0.0 
        do n = 1, nstep_pdd
            day = 1 + (n-1)*dt_pdd
            k1  = frc%k_pdd(n)
            now%t2m = var_today(frc%days(k1-1),frc%days(k1),frc%t2m(i,j,k1-1),frc%t2m(i,j,k1),day)
            call calc_temp_effective(teff,now%t2m-273.15,now%sigma)
            now%PDDs = now%PDDs + teff*real(dt_pdd,prec)
        end do 

        ! Initialize averaging 
        ann = smbpal_point_class()
        if (calc_monthly) mon = smbpal_point_class()

        mnow = 1 
        mday = 0 

        do n = 1, nstep_itm 

            ! Determine t2m, pr, sf and S today 
            day = 1 + (n-1)*dt_itm
            k1  = frc%k_step(n)
            now%t2m = var_today(frc%days(k1-1),frc%days(k1),frc%t2m(i,j,k1-1),frc%t2m(i,j,k1),day)
            now%pr  = var_today(frc%days(k1-1),frc%days(k1),frc%pr(i,j,k1-1), frc%pr(i,j,k1),day)
            now%sf  = var_today(frc%days(k1-1),frc%days(k1),frc%sf(i,j,k1-1), frc%sf(i,j,k1),day)
            
            now%S   = insol_day_eval(frc%insol(n),dble(par%lats(i,j)))

            ! Call mass budget for today [mm/d]
            call calc_snowpack_budget_step(par%itm,dt,par%lats(i,j),z_srf,H_ice,now%S,now%t2m,now%PDDs, &
                                           now%pr,now%sf,now%H_snow,now%alb_s,now%smbi, &
                                           now%smb,now%melt,now%runoff,now%refrz,now%melt_net)

            ! Sum for averages 
            ann = ann + now 

            if (calc_monthly) then 
                mon(mnow) = mon(mnow) + now 

                mday = mday + dt_itm 
                if (mday .eq. ndays_mon) then 
                    mon(mnow) = mon(mnow) / (real(ndays_mon)/dt)
                    mnow = mnow + 1
                    mday = 0 
                end if 
            end if 
    
        end do 

        ! Finalize annual average 
        ann = ann / (real(ndays)/dt)

        ! Convert mass quantities [mm/d] => [mm/a] 
        ann%pr       = ann%pr       *real(ndays)
        ann%sf       = ann%sf       *real(ndays)
        ann%melt     = ann%melt     *real(ndays)
        ann%runoff   = ann%runoff   *real(ndays)
        ann%refrz    = ann%refrz    *real(ndays)
        ann%smb      = ann%smb      *real(ndays)
        ann%smbi     = ann%smbi     *real(ndays)
        ann%melt_net = ann%melt_net *real(ndays)

        ! Calculate surface temp 
        ann%tsrf = calc_temp_surf(ann%t2m,H_ice,ann%melt_net,fac=par%firn_fac)

        return 

    end subroutine smbpal_itm_point

    subroutine smbpal_update_pdd(ann,par,PDDs_ann,z_srf,H_ice,t2m_ann,pr_ann,sf_ann)

        implicit none 

        type(smbpal_state_class), intent(INOUT) :: ann
        type(smbpal_param_class), intent(IN)    :: par 
        real(prec),               intent(IN)    :: PDDs_ann(:,:) 
        real(prec),               intent(IN)    :: z_srf(:,:)
        real(prec),               intent(IN)    :: H_ice(:,:)
        real(prec),               intent(IN)    :: t2m_ann(:,:)
        real(prec),               intent(IN)    :: pr_ann(:,:)
        real(prec),               intent(IN)    :: sf_ann(:,:)

        ! Store known annual values
        ann%PDDs = PDDs_ann 
        ann%t2m  = t2m_ann 
        ann%pr   = pr_ann 
        ann%sf   = sf_ann

        write(*,*) "smbpal_update_pdd"
        write(*,*) "sf:  ", minval(ann%sf), maxval(ann%sf)
        write(*,*) "t2m: ", minval(ann%t2m), maxval(ann%t2m)
        
        ! Get ablation, runoff and refreezing [mm/a]
        call calc_ablation_pdd(ann%melt,ann%runoff,ann%refrz,ann%PDDs,ann%sf, &
                                par%mm_snow,par%mm_ice,par%itm%Pmaxfrac)

        ! Get surface mass balance [mm/a]
        ann%smb  = ann%sf - ann%runoff 
        ann%smbi = ann%smb 

        ! Get melt_net for surface temp calculations [mm/a]
        ann%melt_net = ann%refrz 

        ! Calculate surface temp 
        ann%tsrf = calc_temp_surf(ann%t2m,H_ice,ann%melt_net,fac=par%firn_fac)

        ! Define other missing variables 
        ann%alb_s = 0.0 
        
        return 

    end subroutine smbpal_update_pdd

    subroutine smbpal_end(smbpal)

        implicit none 

        type(smbpal_class) :: smbpal 

        ! Deallocate smbpal state object
        call smbpal_deallocate(smbpal%now)
	
        return 

    end subroutine smbpal_end

! "Fix TSURF calculations!!!" 

    subroutine smbpal_par_load(par,filename,init,group,itm_group)

        type(smbpal_param_class)     :: par
        character(len=*), intent(IN) :: filename 
        logical, optional :: init 
        logical :: init_pars 
        character(len=*),  intent(IN), optional :: group, itm_group

        ! Local variables 
        integer :: file_unit 
        character(len=32) :: nml_group, itm_nml_group

        ! Make sure we know the namelist group for the smbpal block
        if (present(group)) then
            nml_group = trim(group)
        else
            nml_group = "smbpal"         ! Default parameter blcok name
        end if

        ! Make sure we know the namelist group for the itm block
        if (present(itm_group)) then
            itm_nml_group = trim(itm_group)
        else
            itm_nml_group = "itm"         ! Default parameter blcok name
        end if

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE.

        call nml_read(filename,nml_group,"insol_fldr",par%insol_fldr,init=init_pars)
        call nml_read(filename,nml_group,"const_insol",par%const_insol,init=init_pars)
        call nml_read(filename,nml_group,"const_kabp",par%const_kabp,init=init_pars)
        call nml_read(filename,nml_group,"abl_method",par%abl_method,init=init_pars)
        call nml_read(filename,nml_group,"sigma_snow",par%sigma_snow,init=init_pars)
        call nml_read(filename,nml_group,"sigma_land",par%sigma_land,init=init_pars)
        call nml_read(filename,nml_group,"sigma_melt",par%sigma_melt,init=init_pars)
        call nml_read(filename,nml_group,"sf_a",par%sf_a,init=init_pars)
        call nml_read(filename,nml_group,"sf_b",par%sf_b,init=init_pars)
        call nml_read(filename,nml_group,"firn_fac",par%firn_fac,init=init_pars)
        call nml_read(filename,nml_group,"mm_snow",par%mm_snow,init=init_pars)
        call nml_read(filename,nml_group,"mm_ice",par%mm_ice,init=init_pars)

        ! Also load itm parameters
        call itm_par_load(par%itm,filename,init=init,group=itm_nml_group)

        ! Local parameter definitions (identical to object)
        ! character(len=512) :: insol_fldr 
        ! logical    :: const_insol
        ! real(prec) :: const_kabp
        ! character(len=16)  :: abl_method
        ! real(prec)         :: sigma_snow, sigma_melt, sigma_land
        ! real(prec)         :: sf_a, sf_b, firn_fac 
        ! real(prec)         :: mm_snow, mm_ice 

        ! namelist /smbpal/ insol_fldr, const_insol, const_kabp, &
        !     abl_method, sigma_snow, sigma_melt, sigma_land, &
        !     sf_a, sf_b, firn_fac, mm_snow, mm_ice 
                
        ! ! Store initial values in local parameter values 
        ! insol_fldr  = par%insol_fldr
        ! const_insol = par%const_insol
        ! const_kabp  = par%const_kabp
        ! abl_method  = par%abl_method
        ! sigma_snow  = par%sigma_snow 
        ! sigma_melt  = par%sigma_melt 
        ! sigma_land  = par%sigma_land 
        ! sf_a        = par%sf_a 
        ! sf_b        = par%sf_b 
        ! firn_fac    = par%firn_fac 
        ! mm_snow     = par%mm_snow 
        ! mm_ice      = par%mm_ice 

        ! Read parameters from input namelist file
        ! inquire(file=trim(filename),NUMBER=file_unit)
        ! if (file_unit .gt. 0) then 
        !     read(file_unit,nml=smbpal)
        ! else
        !     open(7,file=trim(filename))
        !     read(7,nml=smbpal)
        !     close(7)
        ! end if 

        ! ! Store local parameter values in output object
        ! par%insol_fldr  = insol_fldr 
        ! par%const_insol = const_insol
        ! par%const_kabp  = const_kabp
        ! par%abl_method  = abl_method
        ! par%sigma_snow  = sigma_snow 
        ! par%sigma_melt  = sigma_melt 
        ! par%sigma_land  = sigma_land 
        ! par%sf_a        = sf_a 
        ! par%sf_b        = sf_b 
        ! par%firn_fac    = firn_fac 
        ! par%mm_snow     = mm_snow 
        ! par%mm_ice      = mm_ice 

        ! ! Also load itm parameters
        ! call itm_par_load(par%itm,filename)

        return

    end subroutine smbpal_par_load

   
    ! =======================================================
    !
    ! smb physics (general)
    !
    ! =======================================================

    elemental function calc_temp_surf(tann,H_ice,melt_net,fac) result(ts)
        ! Surface temperature is equal to the annual mean
        ! near-surface temperature + warming due to 
        ! freezing of superimposed ice - cooling due to melt
        implicit none 

        real(prec), intent(IN) :: tann, H_ice, melt_net, fac 
        real(prec) :: ts 

        if (H_ice .gt. 0.0) then
            ! Adjust temp to account for positive melt_net (refreezing) warms firn
            ts = tann + fac * max(0.0, melt_net)
            ! Limit temps to freezing temperature on the ice sheet 
            ts = min(273.15, ts)  
        else
            ts = tann
        end if
        
        return 

    end function calc_temp_surf
    
    elemental function calc_snowfrac(t2m,a,b) result(f)
        ! Return the fraction of snow from total precipitation
        ! expected for a given temperature
        
        implicit none 

        real(prec), intent(IN) :: t2m, a, b 
        real(prec)             :: f 

        f = -0.5*tanh(a*(t2m-b))+0.5 

        return 

    end function calc_snowfrac

    ! =======================================================
    !
    ! smbpal I/O
    !
    ! =======================================================

    subroutine smbpal_write_init(par,filename,z_srf,H_ice)

        implicit none 

        type(smbpal_param_class), intent(IN) :: par 
        character(len=*),         intent(IN) :: filename 
        real(prec), intent(IN), optional :: z_srf(:,:), H_ice(:,:) 

        call nc_create(filename)
        call nc_write_dim(filename,"xc",x=par%x)
        call nc_write_dim(filename,"yc",x=par%y)
        call nc_write_dim(filename,"day",  x=1,nx=360,dx=1)
        call nc_write_dim(filename,"month",x=1,nx=12,dx=1)
        call nc_write_dim(filename,"time",x=0.0,units="kiloyears",unlimited=.TRUE.)
        
        ! Write the 2D latitude field to file
        call nc_write(filename,"lat2D",par%lats,dim1="xc",dim2="yc")

        if (present(z_srf)) call nc_write(filename,"z_srf",z_srf,dim1="xc",dim2="yc")
        if (present(H_ice)) call nc_write(filename,"H_ice",H_ice,dim1="xc",dim2="yc")

        return 

    end subroutine smbpal_write_init

    subroutine smbpal_write(now,filename,time_bp,step,nstep)

        implicit none 

        type(smbpal_state_class), intent(IN) :: now 
        character(len=*),         intent(IN) :: filename 
        real(prec),               intent(IN) :: time_bp  
        character(len=*),         intent(IN) :: step  
        integer, intent(IN), optional        :: nstep 

        ! Local variables 
        real(prec) :: ka_bp 
        integer :: ndat, nx, ny, nt   
        real(prec), allocatable :: time(:) 
        character(len=56) :: step_name 

        ka_bp = time_bp * 1e-3 

        if (trim(step) .ne. "ann" .and. trim(step) .ne. "mon" .and. trim(step) .ne. "day") then 
            write(*,*) "smbpal_write:: error: step should be one of: ann, mon or day."
            stop 
        end if 

        nx = size(now%t2m,1)
        ny = size(now%t2m,2)

        ! Determine timestep to be written 
        nt = nc_size(filename,"time")
        allocate(time(nt))
        call nc_read(filename,"time",time)
        
        if (maxval(time) .lt. ka_bp) then 
            ndat = nt+1 
        else 
            ndat = minloc(abs(time-ka_bp),1)
        end if 

        ! Write the variables
        if (trim(step) .eq. "ann") then 
            ! Write the annual mean with time as 3rd dimension 

            ! Update the timestep 
            call nc_write(filename,"time",ka_bp,dim1="time",start=[ndat],count=[1])

            call nc_write(filename,"t2m",now%t2m,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Near-surface temperature",units="K")
            call nc_write(filename,"S",now%S,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Solar insolation (TOA)",units="W m**-2")
            call nc_write(filename,"pr",now%pr,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Precipitation",units="mm d**-1")
            call nc_write(filename,"sf",now%sf,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Snowfall",units="mm d**-1")
            call nc_write(filename,"PDDs",now%PDDs,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Positive degree days",units="d K")
            call nc_write(filename,"tsrf",now%tsrf,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Ice surface temperature",units="K")

            call nc_write(filename,"H_snow",now%H_snow,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Snowpack thickness",units="mm w.e.")
            call nc_write(filename,"alb_s",now%alb_s,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Surface albedo",units="1")
            call nc_write(filename,"smbi",now%smbi,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Surface mass balance (ice)",units="mm d**-1")
            call nc_write(filename,"smb",now%smb,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Surface mass balance (snow)",units="mm d**-1")
            call nc_write(filename,"melt",now%melt,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Total melt",units="mm d**-1")
            call nc_write(filename,"runoff",now%runoff,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Net runoff",units="mm d**-1")
            call nc_write(filename,"refrz",now%refrz,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,ndat],count=[nx,ny,1],long_name="Refreezing",units="mm d**-1")
            
        else 
            ! Write the step along the year (mon or day)

            if (.not. present(nstep)) then 
                write(*,*) "smbpal_write:: error: nstep must be given for 'day' or 'mon' writing."
                stop 
            end if 

            step_name = "day" 
            if (trim(step) .eq. "mon") step_name = "month" 

            ! Update the timestep 
            call nc_write(filename,"time",ka_bp,dim1="time",start=[1],count=[1])

            call nc_write(filename,"t2m",now%t2m,dim1="xc",dim2="yc",dim3=trim(step_name), &
                          start=[1,1,nstep],count=[nx,ny,1])
            call nc_write(filename,"S",now%S,dim1="xc",dim2="yc",dim3=trim(step_name), &
                          start=[1,1,nstep],count=[nx,ny,1])

        end if 

        return 

    end subroutine smbpal_write

    subroutine smbpal_restart_write(smb,filename,time)
        ! Write the prognostic snowpack state carried across timesteps
        ! (H_snow, alb_s). With the ITM ablation method the snowpack is
        ! prognostic, so it must be restored to continue a run seamlessly.

        implicit none

        type(smbpal_class), intent(IN) :: smb
        character(len=*),   intent(IN) :: filename
        real(prec),         intent(IN) :: time

        ! Local variables
        integer :: ncid, nx, ny

        nx = size(smb%now%H_snow,1)
        ny = size(smb%now%H_snow,2)

        call nc_create(filename)
        call nc_write_dim(filename,"xc",  x=smb%par%x)
        call nc_write_dim(filename,"yc",  x=smb%par%y)
        call nc_write_dim(filename,"time",x=time,dx=1.0_prec,nx=1,units="year",unlimited=.TRUE.)

        call nc_open(filename,ncid,writable=.TRUE.)
        call nc_write(filename,"H_snow",smb%now%H_snow,dim1="xc",dim2="yc",dim3="time", &
                      start=[1,1,1],count=[nx,ny,1],ncid=ncid)
        call nc_write(filename,"alb_s", smb%now%alb_s, dim1="xc",dim2="yc",dim3="time", &
                      start=[1,1,1],count=[nx,ny,1],ncid=ncid)
        call nc_close(ncid)

        return

    end subroutine smbpal_restart_write

    subroutine smbpal_restart_read(smb,filename)
        ! Restore the prognostic snowpack state (H_snow, alb_s).

        implicit none

        type(smbpal_class), intent(INOUT) :: smb
        character(len=*),   intent(IN)    :: filename

        ! Local variables
        integer :: nx, ny

        nx = size(smb%now%H_snow,1)
        ny = size(smb%now%H_snow,2)

        call nc_read(filename,"H_snow",smb%now%H_snow,start=[1,1,1],count=[nx,ny,1])
        call nc_read(filename,"alb_s", smb%now%alb_s, start=[1,1,1],count=[nx,ny,1])

        write(*,*) "smbpal_restart_read:: read "//trim(filename)

        return

    end subroutine smbpal_restart_read

    ! =======================================================
    !
    ! smbpal memory / data management
    !
    ! =======================================================

    subroutine smbpal_allocate(now,nx,ny)

        implicit none 

        type(smbpal_state_class) :: now 
        integer :: nx, ny 

        ! Make object is deallocated
        call smbpal_deallocate(now)

        ! Allocate variables
        allocate(now%t2m(nx,ny))
        allocate(now%pr(nx,ny))
        allocate(now%sf(nx,ny))
        allocate(now%S(nx,ny))
        allocate(now%sigma(nx,ny))
        allocate(now%PDDs(nx,ny))
        allocate(now%tsrf(nx,ny))
        allocate(now%H_snow(nx,ny))
        allocate(now%alb_s(nx,ny))
        allocate(now%smbi(nx,ny))
        allocate(now%smb(nx,ny))
        allocate(now%melt(nx,ny))
        allocate(now%runoff(nx,ny))
        allocate(now%refrz(nx,ny))

        allocate(now%melt_net(nx,ny))

        ! Define every field: some are only set by one ablation method (e.g.
        ! alb_s by itm), but all are written to output and restart files.
        now%t2m      = 0.0
        now%pr       = 0.0
        now%sf       = 0.0
        now%S        = 0.0
        now%sigma    = 0.0
        now%PDDs     = 0.0
        now%tsrf     = 0.0
        now%H_snow   = 0.0
        now%alb_s    = 0.0
        now%smbi     = 0.0
        now%smb      = 0.0
        now%melt     = 0.0
        now%runoff   = 0.0
        now%refrz    = 0.0
        now%melt_net = 0.0

        return

    end subroutine smbpal_allocate

    subroutine smbpal_deallocate(now)

        implicit none 

        type(smbpal_state_class) :: now 

        ! Allocate state objects
        if (allocated(now%t2m))      deallocate(now%t2m)
        if (allocated(now%pr))       deallocate(now%pr)
        if (allocated(now%sf))       deallocate(now%sf)
        if (allocated(now%S))        deallocate(now%S)
        if (allocated(now%sigma))    deallocate(now%sigma)
        if (allocated(now%PDDs))     deallocate(now%PDDs)
        if (allocated(now%tsrf))     deallocate(now%tsrf)
        if (allocated(now%H_snow))   deallocate(now%H_snow)
        if (allocated(now%alb_s))    deallocate(now%alb_s)
        if (allocated(now%smbi))     deallocate(now%smbi)
        if (allocated(now%smb))      deallocate(now%smb)
        if (allocated(now%melt))     deallocate(now%melt)
        if (allocated(now%runoff))   deallocate(now%runoff)
        if (allocated(now%refrz))    deallocate(now%refrz)
        
        if (allocated(now%melt_net))   deallocate(now%melt_net)
        
        return

    end subroutine smbpal_deallocate

    function smbpal_point_get(st,i,j) result(p)
        ! State at grid point (i,j)
        implicit none 

        type(smbpal_state_class), intent(IN) :: st 
        integer, intent(IN) :: i, j 
        type(smbpal_point_class) :: p 

        p%t2m      = st%t2m(i,j)
        p%pr       = st%pr(i,j)
        p%sf       = st%sf(i,j)
        p%S        = st%S(i,j)
        p%sigma    = st%sigma(i,j)
        p%PDDs     = st%PDDs(i,j)
        p%tsrf     = st%tsrf(i,j)
        p%H_snow   = st%H_snow(i,j)
        p%alb_s    = st%alb_s(i,j)
        p%smbi     = st%smbi(i,j)
        p%smb      = st%smb(i,j)
        p%melt     = st%melt(i,j)
        p%runoff   = st%runoff(i,j)
        p%refrz    = st%refrz(i,j)
        p%melt_net = st%melt_net(i,j)

        return

    end function smbpal_point_get

    subroutine smbpal_point_set(st,i,j,p)
        ! Store the state of grid point (i,j)
        implicit none 

        type(smbpal_state_class), intent(INOUT) :: st 
        integer, intent(IN) :: i, j 
        type(smbpal_point_class), intent(IN) :: p 

        st%t2m(i,j)      = p%t2m
        st%pr(i,j)       = p%pr
        st%sf(i,j)       = p%sf
        st%S(i,j)        = p%S
        st%sigma(i,j)    = p%sigma
        st%PDDs(i,j)     = p%PDDs
        st%tsrf(i,j)     = p%tsrf
        st%H_snow(i,j)   = p%H_snow
        st%alb_s(i,j)    = p%alb_s
        st%smbi(i,j)     = p%smbi
        st%smb(i,j)      = p%smb
        st%melt(i,j)     = p%melt
        st%runoff(i,j)   = p%runoff
        st%refrz(i,j)    = p%refrz
        st%melt_net(i,j) = p%melt_net

        return

    end subroutine smbpal_point_set

    elemental function smbpal_point_add(a,b) result(c)
        ! Field-wise sum of two point states (for time averages)
        implicit none 

        type(smbpal_point_class), intent(IN) :: a, b 
        type(smbpal_point_class) :: c 

        c%t2m      = a%t2m      + b%t2m
        c%pr       = a%pr       + b%pr
        c%sf       = a%sf       + b%sf
        c%S        = a%S        + b%S
        c%sigma    = a%sigma    + b%sigma
        c%PDDs     = a%PDDs     + b%PDDs
        c%tsrf     = a%tsrf     + b%tsrf
        c%H_snow   = a%H_snow   + b%H_snow
        c%alb_s    = a%alb_s    + b%alb_s
        c%smbi     = a%smbi     + b%smbi
        c%smb      = a%smb      + b%smb
        c%melt     = a%melt     + b%melt
        c%runoff   = a%runoff   + b%runoff
        c%refrz    = a%refrz    + b%refrz
        c%melt_net = a%melt_net + b%melt_net

        return

    end function smbpal_point_add

    elemental function smbpal_point_div(a,nt) result(c)
        ! Field-wise division of a point state (for time averages)
        implicit none 

        type(smbpal_point_class), intent(IN) :: a 
        real(prec), intent(IN) :: nt 
        type(smbpal_point_class) :: c 

        c%t2m      = a%t2m      / nt
        c%pr       = a%pr       / nt
        c%sf       = a%sf       / nt
        c%S        = a%S        / nt
        c%sigma    = a%sigma    / nt
        c%PDDs     = a%PDDs     / nt
        c%tsrf     = a%tsrf     / nt
        c%H_snow   = a%H_snow   / nt
        c%alb_s    = a%alb_s    / nt
        c%smbi     = a%smbi     / nt
        c%smb      = a%smb      / nt
        c%melt     = a%melt     / nt
        c%runoff   = a%runoff   / nt
        c%refrz    = a%refrz    / nt
        c%melt_net = a%melt_net / nt

        return

    end function smbpal_point_div

    function idx_today(days,day) result(idx)

        implicit none 

        integer :: days(:), day
        integer :: idx 

        ! Determine the index of today 
        do idx = 2, size(days)
            if (days(idx) .ge. day) exit
        end do 

        return 

    end function idx_today 
   
    elemental function var_today(x0,x1,y0,y1,x) result(y)
        ! Interpolate y0 and y1 to y (can be fields) assuming that 
        ! x lies within x0 and x1
        implicit none 
        integer, intent(IN) :: x0, x1, x
        real(prec), intent(IN) :: y0, y1
        real(prec) :: y 
        real(prec) :: alpha 

        alpha = dble(x - x0) / dble(x1 - x0)
        y     = y0 + alpha*(y1-y0)

        return 

    end function var_today

end module smbpal


