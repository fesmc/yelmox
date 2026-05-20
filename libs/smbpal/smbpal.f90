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
        real(prec) :: rho_sw
        real(prec) :: rho_ice

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
    
    private
    public :: smbpal_class
    public :: smbpal_init 
    public :: smbpal_update_2temp, smbpal_update_monthly 
    public :: smbpal_update_monthly_equil
    public :: smbpal_end 
    public :: smbpal_write_init, smbpal_write

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
                                   file_out,file_out_mon,file_out_day,write_init,calc_mon,write_now)
        ! Generate climate using two points in year (Tsum,Tann)

        implicit none 
        
        type(smbpal_class), intent(INOUT) :: smb
        real(prec), intent(IN) :: t2m_ann(:,:), t2m_sum(:,:)
        real(prec), intent(IN) ::  pr_ann(:,:), z_srf(:,:), H_ice(:,:)
        real(prec), intent(IN) :: time_bp       ! years BP 
        real(prec), intent(IN), optional :: sf_ann(:,:)
        character(len=*), intent(IN), optional :: file_out      ! Annual output
        character(len=*), intent(IN), optional :: file_out_mon  ! Monthly output
        character(len=*), intent(IN), optional :: file_out_day  ! Daily output 
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
                        file_out,file_out_mon,file_out_day,write_init,calc_mon,write_now)
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

        ! Loop over equilibration years to update snowpack thickness 
        do n = 1, int(time_equil) 
            call smbpal_update_monthly(smb,t2m,pr,z_srf,H_ice,time_bp,sf)
        end do 


        return 

    end subroutine smbpal_update_monthly_equil


    subroutine smbpal_update_monthly(smb,t2m,pr,z_srf,H_ice,time_bp,sf, &
                        file_out,file_out_mon,file_out_day,write_init,calc_mon,write_now)
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
        character(len=*),   intent(IN), optional :: file_out_day    ! Daily output filename 
        logical,            intent(IN), optional :: write_init      ! Flag for whether to initialize writing of output file
        logical,            intent(IN), optional :: calc_mon        ! Flag for whether to calculate monthly averages (itm only)
        logical,            intent(IN), optional :: write_now       ! Flag for whether to write the current time to file

        ! Local variables
        logical :: init_now, write_out_now
        integer :: ndays_daily, k, day, k1   
        integer, allocatable :: daily(:)
        real(prec), allocatable :: t2m_daily(:,:,:), pr_daily(:,:,:)
        real(prec), allocatable :: sf_daily(:,:,:)
        double precision, allocatable :: tmp(:,:,:)
        
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
            ndays_daily = 37
            allocate(daily(ndays_daily))
            allocate(t2m_daily(size(t2m,1),size(t2m,2),ndays_daily))
            allocate(pr_daily(size(t2m,1),size(t2m,2),ndays_daily))
            allocate(sf_daily(size(t2m,1),size(t2m,2),ndays_daily))
            allocate(tmp(size(t2m,1),size(t2m,2),ndays_daily))
            
            ! Define daily days
            do k = 1, ndays_daily-1 
                daily(k) = 1 + (k-1)*(ndays / (ndays_daily-1))
            end do 
            daily(ndays_daily) = ndays 

            ! Generate daily climate from monthly input 
            call convert_monthly_daily_3D(dble(t2m),tmp,days=daily)
            t2m_daily = tmp 
            call convert_monthly_daily_3D(dble(pr),tmp,days=daily)
            pr_daily = tmp 

            if (present(sf)) then 
                call convert_monthly_daily_3D(dble(sf),tmp,days=daily)
                sf_daily = tmp 
                where(sf_daily .lt. 0.0) sf_daily = 0.0 
            else 
                do k = 1, ndays_daily 
                    sf_daily(:,:,k) = pr_daily(:,:,k) * calc_snowfrac(t2m_daily(:,:,k),smb%par%sf_a,smb%par%sf_b)
                end do 
            end if 

                
            ! Call daily subroutine 
            call smbpal_update_itm(smb,daily,t2m_daily,pr_daily,sf_daily,z_srf,H_ice,time_bp, &
                                   file_out_mon,file_out_day,write_init,calc_mon,write_now)
        
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

    subroutine smbpal_update_itm(smb,days,t2m,pr,sf,z_srf,H_ice,time_bp, &
                                   file_out_mon,file_out_day,write_init,calc_mon,write_now)
        ! Generate smb using daily input of climate [nx,ny,nday]

        implicit none 
        
        type(smbpal_class), intent(INOUT) :: smb
        integer, intent(IN) :: days(:)
        real(prec), intent(IN) :: t2m(:,:,:), pr(:,:,:), sf(:,:,:)
        real(prec), intent(IN) ::  z_srf(:,:), H_ice(:,:)
        real(prec), intent(IN) :: time_bp       ! years BP
        character(len=*), intent(IN), optional :: file_out_mon  ! Monthly output
        character(len=*), intent(IN), optional :: file_out_day  ! Daily output 
        logical, intent(IN), optional :: write_init, calc_mon, write_now

        ! Local variables
        logical :: init_now, calc_monthly, write_out_now   
        integer, parameter :: ndays = 360       ! 360-day year
        integer, parameter :: ndays_mon = 30    ! 30 days per month  
        integer :: day, m, nx, ny, mnow, mday  
        integer :: k1 
        real(prec) :: dt    ! [days]
        real(8) :: insol_time
        
        type(smbpal_param_class) :: par
        type(smbpal_state_class) :: now

        real(prec), allocatable :: tmp(:,:) 

        allocate(tmp(size(t2m,1),size(t2m,2)))

        ! Determine whether this is first time running (for output)
        init_now = .FALSE. 
        if (present(write_init)) init_now = write_init 

        calc_monthly = .FALSE. 
        if (present(calc_mon))     calc_monthly = calc_mon 
        if (present(file_out_mon)) calc_monthly = .TRUE. 

        write_out_now = .FALSE. 
        if (present(write_now)) write_out_now = write_now 
        
        ! Determine year to use for insolation calcs
        insol_time = time_bp
        if (smb%par%const_insol) insol_time = smb%par%const_kabp*1e3
        
        ! Set sigma to snow sigma everywhere for pdd calcs
        smb%now%sigma = smb%par%sigma_snow
        
        ! Fill in local versions for easier access 
        par = smb%par 
        now = smb%now 

        ! First calculate PDDs for the whole year (input to itm)
        now%PDDs = 0.0 
        do day = 1, ndays, 10
            k1 = idx_today(days,day)
            now%t2m = var_today(days(k1-1),days(k1),t2m(:,:,k1-1),t2m(:,:,k1),day)
            call calc_temp_effective(tmp,now%t2m-273.15,now%sigma)
            now%PDDs = now%PDDs + tmp*10.0
        end do 

        ! Initialize averaging 
        call smbpal_average(smb%ann,now,step="init")

        if (calc_monthly) then 
            do m = 1, 12 
                call smbpal_average(smb%mon(m),now,step="init")
            end do
        end if 

        mnow = 1 
        mday = 0 

        ! Initialize daily output file if needed
        if (write_out_now .and. present(file_out_day)) then
            if (init_now) call smbpal_write_init(par,file_out_day,z_srf,H_ice)
        end if 

        dt = 2.0 

        do day = 1, ndays, int(dt)

            ! Determine t2m, pr, sf and S today 
            k1 = idx_today(days,day)
            now%t2m = var_today(days(k1-1),days(k1),t2m(:,:,k1-1),t2m(:,:,k1),day)
            now%pr  = var_today(days(k1-1),days(k1),pr(:,:,k1-1), pr(:,:,k1),day)
            now%sf  = var_today(days(k1-1),days(k1),sf(:,:,k1-1), sf(:,:,k1),day)
            
            now%S   = calc_insol_day(day,dble(par%lats),insol_time,fldr=par%insol_fldr)

            ! Call mass budget for today [mm/d]
            call calc_snowpack_budget_step(par%itm,dt,par%lats,z_srf,H_ice,now%S,now%t2m,now%PDDs, &
                                           now%pr,now%sf,now%H_snow,now%alb_s,now%smbi, &
                                           now%smb,now%melt,now%runoff,now%refrz,now%melt_net)
        

            ! Get averages 
            call smbpal_average(smb%ann,now,step="step")

            if (calc_monthly) then 
                call smbpal_average(smb%mon(mnow),now,step="step")

                mday = mday + int(dt) 
                if (mday .eq. ndays_mon) then 
                    call smbpal_average(smb%mon(mnow),now,step="end",nt=real(ndays_mon)/dt)
                    mnow = mnow + 1
                    mday = 0 
                end if 
            end if 

            if (write_out_now .and. present(file_out_day)) then 
                ! Write daily output for this year 
                call smbpal_write(now,file_out_day,time_bp=time_bp,step="day",nstep=day)
            end if 
    
        end do 

        ! Finalize annual average 
        call smbpal_average(smb%ann,now,step="end",nt=real(ndays)/dt)

        ! Convert mass quantities [mm/d] => [mm/a] 
        smb%ann%pr       = smb%ann%pr       *real(ndays)
        smb%ann%sf       = smb%ann%sf       *real(ndays)
        smb%ann%melt     = smb%ann%melt     *real(ndays)
        smb%ann%runoff   = smb%ann%runoff   *real(ndays)
        smb%ann%refrz    = smb%ann%refrz    *real(ndays)
        smb%ann%smb      = smb%ann%smb      *real(ndays)
        smb%ann%smbi     = smb%ann%smbi     *real(ndays)
        smb%ann%melt_net = smb%ann%melt_net *real(ndays)

        ! Calculate surface temp 
        smb%ann%tsrf = calc_temp_surf(smb%ann%t2m,H_ice,smb%ann%melt_net,fac=par%firn_fac)

        ! Repopulate global now variable (in case it is needed)
        smb%now = now 

        ! Monthly I/O 
        if (write_out_now .and. calc_monthly .and. present(file_out_mon)) then
            if (init_now) call smbpal_write_init(par,file_out_mon,z_srf,H_ice)
            do m = 1, 12
                call smbpal_write(smb%mon(m),file_out_mon,time_bp=time_bp,step="mon",nstep=m)
            end do 

        end if 

        return 

    end subroutine smbpal_update_itm

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

    subroutine smbpal_average(ave,now,step,nt)
        implicit none 

        type(smbpal_state_class), intent(INOUT) :: ave
        type(smbpal_state_class), intent(IN)    :: now 
        character(len=*)  :: step
        real(prec), optional :: nt 
        
        call field_average(ave%t2m,    now%t2m,    step,nt)
        call field_average(ave%pr,     now%pr,     step,nt)
        call field_average(ave%sf,     now%sf,     step,nt)
        call field_average(ave%S,      now%S,      step,nt)
        
        call field_average(ave%H_snow, now%H_snow, step,nt)
        call field_average(ave%alb_s,  now%alb_s,  step,nt)
        call field_average(ave%smbi,   now%smbi,   step,nt)
        call field_average(ave%smb,    now%smb,    step,nt)
        call field_average(ave%melt,   now%melt,   step,nt)
        call field_average(ave%runoff, now%runoff, step,nt)
        call field_average(ave%refrz,  now%refrz,  step,nt)
        
        call field_average(ave%melt_net,now%melt_net,step,nt)
        
        ! Annual values, averaged for completeness 
        call field_average(ave%sigma,  now%sigma,  step,nt)
        call field_average(ave%PDDs,   now%PDDs,   step,nt)
        call field_average(ave%tsrf,   now%tsrf,   step,nt)
        
        return

    end subroutine smbpal_average

    subroutine field_average(ave,now,step,nt)
        ! Generic routine to average a field through time 

        implicit none 
        real(prec), intent(INOUT)    :: ave(:,:)
        real(prec), intent(IN)       :: now(:,:)
        character(len=*), intent(IN) :: step
        real(prec), intent(IN), optional :: nt 

        if (trim(step) .eq. "init") then
            ! Initialize field to zero  
            ave = 0.0 
        else if (trim(step) .eq. "step") then 
            ! Sum intermediate steps
            ave = ave + now 
        else if (trim(step) .eq. "end") then
            if (.not.  present(nt)) then 
                write(*,*) "Averaging step total not provided."
                stop 
            end if 
            ! Divide by total steps
            ave = ave / nt 
        else
            write(*,*) "Step not recognized: ",trim(step)
            stop 
        end if 

        return 

    end subroutine field_average 

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


