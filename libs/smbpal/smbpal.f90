 module smbpal

    use nml
    use smbpal_precision
    use insolation, only : insol_class, insol_init, calc_insol_days, insol_wp => wp
    use orbit_forcing, only : orbit_ctl_class, orbit_pars_class, orbit_ctl_load, &
                              orbit_pars_now, orbit_pars_equal
    use interp_time 
    use ncio
    use smb_pdd 
    use smb_itm 

    implicit none 

    integer, parameter :: ndays = 360   ! 360-day year

    ! itm: daily forcing interpolated from the monthly fields at ndays_daily
    ! days, snowpack budget every dt_itm days, annual PDDs from every dt_pdd days
    integer, parameter :: ndays_daily = 37
    integer, parameter :: dt_itm      = 2
    integer, parameter :: nstep_itm   = ndays/dt_itm
    integer, parameter :: dt_pdd      = 10
    integer, parameter :: nstep_pdd   = ndays/dt_pdd
        
    type smbpal_param_class
        type(itm_par_class) :: itm
        type(orbit_ctl_class) :: orbit        ! [orbit<sfx>]: the orbit of the insolation
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

    type smbpal_insol_class
        ! Daily TOA insolation of each itm step on the grid (itm only),
        ! recomputed only when the orbital elements change.
        type(insol_class)       :: is                   ! insol state on the grid latitudes
        type(orbit_pars_class)  :: orbit                ! the elements of S
        logical                 :: set = .FALSE.        ! S holds the insolation of orbit
        real(prec), allocatable :: S(:,:,:)             ! [W/m2] [nstep_itm,nx,ny]
    end type

    type smbpal_class
        type(smbpal_param_class) :: par 
        type(smbpal_state_class) :: now, ann
        type(smbpal_insol_class) :: ins
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
    public :: smbpal_update_monthly
    public :: smbpal_update_monthly_equil
    public :: smbpal_restart_write, smbpal_restart_read

contains 

    subroutine smbpal_init(smb,filename,x,y,lats,group,itm_group,orbit_group)

        implicit none 

        type(smbpal_class) :: smb
        character(len=*), intent(IN)  :: filename  ! Parameter file 
        real(prec) :: x(:), y(:), lats(:,:)
        character(len=*),  intent(IN), optional :: group, itm_group, orbit_group

        ! Local variables
        integer :: nx, ny
        character(len=32) :: nml_group, itm_nml_group, orbit_nml_group

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

        ! Make sure we know the namelist group of the orbit
        if (present(orbit_group)) then
            orbit_nml_group = trim(orbit_group)
        else
            orbit_nml_group = "orbit"     ! Default parameter block name
        end if

        nx = size(x,1)
        ny = size(y,1)

        ! Load smbpal parameters
        call smbpal_par_load(smb%par,filename,init=.TRUE.,group=nml_group,itm_group=itm_nml_group)
        call orbit_ctl_load(smb%par%orbit,filename,orbit_nml_group)

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

        ! Initialize the state variables 
        smb%now%H_snow = smb%par%itm%H_snow_max 

        ! Insolation state on the grid latitudes (itm only)
        if (trim(smb%par%abl_method) .eq. "itm") then 
            call insol_init(smb%ins%is,dble(lats),fldr=trim(smb%par%orbit%fldr),day_year=ndays)
            allocate(smb%ins%S(nstep_itm,nx,ny))
            smb%ins%S   = 0.0 
            smb%ins%set = .FALSE. 
        end if 

        return 

    end subroutine smbpal_init

    subroutine smbpal_update_monthly_equil(smb,t2m,pr,z_srf,H_ice,time_bp,time_equil)
        ! Generate climate using monthly input data [nx,ny,nmon]
        
        implicit none 
        
        type(smbpal_class), intent(INOUT) :: smb
        real(prec), intent(IN) :: t2m(:,:,:), pr(:,:,:)
        real(prec), intent(IN) ::  z_srf(:,:), H_ice(:,:)
        real(dp),   intent(IN) :: time_bp       ! [yr] model time (relative to 1950, past < 0), for the orbit 
        real(prec), intent(IN) :: time_equil    ! years to equilibrate

        ! Local variables
        integer :: n 
        type(smbpal_itm_forcing_class) :: frc 

        ! Loop over equilibration years to update snowpack thickness 
        if (trim(smb%par%abl_method) .eq. "itm") then 
            ! The daily forcing is the same every year: prepare it once
            call smbpal_insol_update(smb%ins,orbit_pars_now(smb%par%orbit,smb%ins%is,time_bp))
            call smbpal_itm_forcing_init(frc,smb%par,t2m,pr)
            do n = 1, int(time_equil)
                call smbpal_update_itm(smb,frc,z_srf,H_ice)
            end do
        else 
            do n = 1, int(time_equil) 
                call smbpal_update_monthly(smb,t2m,pr,z_srf,H_ice,time_bp)
            end do 
        end if 

        return 

    end subroutine smbpal_update_monthly_equil


    subroutine smbpal_update_monthly(smb,t2m,pr,z_srf,H_ice,time_bp)
        ! Generate climate using monthly input data [nx,ny,nmon]

        implicit none

        type(smbpal_class), intent(INOUT) :: smb
        real(prec),         intent(IN) :: t2m(:,:,:)                ! [K] Monthly temperature fields
        real(prec),         intent(IN) :: pr(:,:,:)                 ! [mm we/d] Monthly precipitation rate fields
        real(prec),         intent(IN) :: z_srf(:,:)                ! [m] Surface elevation
        real(prec),         intent(IN) :: H_ice(:,:)                ! [m] Ice thickness
        real(dp),           intent(IN) :: time_bp                   ! [yr] Model time (relative to 1950, past < 0), for the orbit

        ! Local variables
        integer :: k
        type(smbpal_itm_forcing_class) :: frc

        real(prec), allocatable :: tmp4(:,:)
        real(prec), allocatable :: t2m_ann(:,:), pr_ann(:,:), sf_ann(:,:)
        real(prec), allocatable :: PDDs_ann(:,:)

        allocate(tmp4(size(t2m,1),size(t2m,2)))

        if (trim(smb%par%abl_method) .eq. "itm") then

            ! Generate daily climate from monthly input and run the year
            call smbpal_insol_update(smb%ins,orbit_pars_now(smb%par%orbit,smb%ins%is,time_bp))
            call smbpal_itm_forcing_init(frc,smb%par,t2m,pr)
            call smbpal_update_itm(smb,frc,z_srf,H_ice)

        else
            ! PDD method 

            allocate(t2m_ann(size(t2m,1),size(t2m,2)))
            allocate(pr_ann(size(t2m,1),size(t2m,2)))
            allocate(sf_ann(size(t2m,1),size(t2m,2)))
            allocate(PDDs_ann(size(t2m,1),size(t2m,2)))

            t2m_ann = sum(t2m,dim=3) / 12.0 
            pr_ann  = sum(pr, dim=3) / 12.0 *real(ndays,prec)       ! [mm we/d] => [mm we/a]

            sf_ann  = sum(calc_snowfall(smb%par,t2m,pr),dim=3) / 12.0 *real(ndays,prec)   ! [mm we/d] => [mm we/a]

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

        return

    end subroutine smbpal_update_monthly

    subroutine smbpal_insol_update(ins,orbit)
        ! Daily TOA insolation of each itm step on the grid, only when the
        ! orbital elements have changed since the last call.

        implicit none 

        type(smbpal_insol_class), intent(INOUT) :: ins 
        type(orbit_pars_class),   intent(IN)    :: orbit 

        ! Local variables
        integer :: n 
        integer :: days(nstep_itm)
        real(insol_wp), allocatable :: S(:,:,:)

        if (ins%set) then 
            if (orbit_pars_equal(orbit,ins%orbit)) return 
        end if 

        do n = 1, nstep_itm 
            days(n) = 1 + (n-1)*dt_itm
        end do 

        allocate(S(size(ins%S,2),size(ins%S,3),nstep_itm))
        call calc_insol_days(ins%is,S,days,orbit)
        do n = 1, nstep_itm 
            ins%S(n,:,:) = real(S(:,:,n),prec)
        end do 

        ins%orbit = orbit 
        ins%set   = .TRUE. 

        return 

    end subroutine smbpal_insol_update

    subroutine smbpal_itm_forcing_init(frc,par,t2m,pr)
        ! Daily forcing of one itm year: the monthly fields interpolated to
        ! the forcing days. It depends only on the climate, so it is shared
        ! by all grid points and years.

        implicit none 

        type(smbpal_itm_forcing_class), intent(OUT) :: frc 
        type(smbpal_param_class),       intent(IN)  :: par 
        real(prec), intent(IN) :: t2m(:,:,:)                ! [K] Monthly temperature fields
        real(prec), intent(IN) :: pr(:,:,:)                 ! [mm we/d] Monthly precipitation rate fields 

        ! Local variables
        integer :: nx, ny, k, n 
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

        ! Snowfall from the daily temperature and precipitation
        frc%sf = calc_snowfall(par,frc%t2m,frc%pr)

        return 

    end subroutine smbpal_itm_forcing_init

    subroutine smbpal_update_itm(smb,frc,z_srf,H_ice)
        ! One itm year from the daily forcing frc. Grid points are
        ! independent, so the year is integrated point by point in parallel.

        implicit none

        type(smbpal_class),             intent(INOUT) :: smb
        type(smbpal_itm_forcing_class), intent(IN)    :: frc
        real(prec),                     intent(IN)    :: z_srf(:,:), H_ice(:,:)

        ! Local variables
        integer :: i, j, nx, ny
        type(smbpal_point_class) :: now, ann

        nx = size(z_srf,1)
        ny = size(z_srf,2)

        !$omp parallel do collapse(2) private(i,j,now,ann)
        do j = 1, ny
        do i = 1, nx

            now = smbpal_point_get(smb%now,i,j)

            call smbpal_itm_point(now,ann,smb%par,frc,smb%ins%S(:,i,j),i,j,z_srf(i,j),H_ice(i,j))

            call smbpal_point_set(smb%now,i,j,now)
            call smbpal_point_set(smb%ann,i,j,ann)

        end do
        end do 
        !$omp end parallel do

        return 

    end subroutine smbpal_update_itm

    subroutine smbpal_itm_point(now,ann,par,frc,S,i,j,z_srf,H_ice)
        ! One itm year at grid point (i,j): the annual PDDs, then the
        ! snowpack budget every dt_itm days, averaged over the year.

        implicit none

        type(smbpal_point_class),       intent(INOUT) :: now        ! State, carried between years
        type(smbpal_point_class),       intent(INOUT) :: ann        ! Annual mean
        type(smbpal_param_class),       intent(IN)    :: par
        type(smbpal_itm_forcing_class), intent(IN)    :: frc
        real(prec),                     intent(IN)    :: S(nstep_itm)   ! [W/m2] TOA insolation of each itm step
        integer,                        intent(IN)    :: i, j
        real(prec),                     intent(IN)    :: z_srf, H_ice

        ! Local variables
        integer :: n, day, k1
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

        do n = 1, nstep_itm

            ! Determine t2m, pr, sf and S today 
            day = 1 + (n-1)*dt_itm
            k1  = frc%k_step(n)
            now%t2m = var_today(frc%days(k1-1),frc%days(k1),frc%t2m(i,j,k1-1),frc%t2m(i,j,k1),day)
            now%pr  = var_today(frc%days(k1-1),frc%days(k1),frc%pr(i,j,k1-1), frc%pr(i,j,k1),day)
            now%sf  = var_today(frc%days(k1-1),frc%days(k1),frc%sf(i,j,k1-1), frc%sf(i,j,k1),day)
            
            now%S   = S(n)

            ! Call mass budget for today [mm/d]
            call calc_snowpack_budget_step(par%itm,dt,par%lats(i,j),z_srf,H_ice,now%S,now%t2m,now%PDDs, &
                                           now%pr,now%sf,now%H_snow,now%alb_s,now%smbi, &
                                           now%smb,now%melt,now%runoff,now%refrz,now%melt_net)

            ! Sum for averages
            ann = ann + now

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

    subroutine smbpal_par_load(par,filename,init,group,itm_group)

        type(smbpal_param_class)     :: par
        character(len=*), intent(IN) :: filename
        logical, optional :: init
        logical :: init_pars
        character(len=*),  intent(IN), optional :: group, itm_group

        ! Local variables
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
        if (present(init)) init_pars = init

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
        call itm_par_load(par%itm,filename,init=init_pars,group=itm_nml_group)

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
    
    elemental function calc_snowfall(par,t2m,pr) result(sf)
        ! Snowfall from precipitation and temperature: the one place
        ! smbpal splits precipitation into snow and rain (itm on the
        ! daily fields, pdd on the monthly fields).

        implicit none 

        type(smbpal_param_class), intent(IN) :: par 
        real(prec),               intent(IN) :: t2m     ! [K]
        real(prec),               intent(IN) :: pr      ! [mm we/d]
        real(prec)                           :: sf      ! [mm we/d]

        sf = pr * calc_snowfrac(t2m,par%sf_a,par%sf_b)

        return 

    end function calc_snowfall

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

    subroutine smbpal_restart_write(smb,filename,time)
        ! Write the prognostic snowpack state carried across timesteps
        ! (H_snow, alb_s). With the ITM ablation method the snowpack is
        ! prognostic, so it must be restored to continue a run seamlessly.
        ! Also the annual fields of the last update (smb%ann, as ann_<name>):
        ! they are the output the host receives at its next step and cannot
        ! be rebuilt from the state without stepping it another year.

        implicit none

        type(smbpal_class), intent(IN) :: smb
        character(len=*),   intent(IN) :: filename
        real(dp),           intent(IN) :: time

        ! Local variables
        integer :: ncid, nx, ny

        nx = size(smb%now%H_snow,1)
        ny = size(smb%now%H_snow,2)

        call nc_create(filename)
        call nc_write_dim(filename,"xc",  x=smb%par%x)
        call nc_write_dim(filename,"yc",  x=smb%par%y)
        call nc_write_dim(filename,"time",x=time,dx=1.0_dp,nx=1,units="year",unlimited=.TRUE.)

        call nc_open(filename,ncid,writable=.TRUE.)
        call nc_write(filename,"H_snow",smb%now%H_snow,dim1="xc",dim2="yc",dim3="time", &
                      start=[1,1,1],count=[nx,ny,1],ncid=ncid)
        call nc_write(filename,"alb_s", smb%now%alb_s, dim1="xc",dim2="yc",dim3="time", &
                      start=[1,1,1],count=[nx,ny,1],ncid=ncid)

        call write_ann("t2m",      smb%ann%t2m)
        call write_ann("pr",       smb%ann%pr)
        call write_ann("sf",       smb%ann%sf)
        call write_ann("S",        smb%ann%S)
        call write_ann("sigma",    smb%ann%sigma)
        call write_ann("PDDs",     smb%ann%PDDs)
        call write_ann("tsrf",     smb%ann%tsrf)
        call write_ann("H_snow",   smb%ann%H_snow)
        call write_ann("alb_s",    smb%ann%alb_s)
        call write_ann("smbi",     smb%ann%smbi)
        call write_ann("smb",      smb%ann%smb)
        call write_ann("melt",     smb%ann%melt)
        call write_ann("runoff",   smb%ann%runoff)
        call write_ann("refrz",    smb%ann%refrz)
        call write_ann("melt_net", smb%ann%melt_net)
        call nc_close(ncid)

        return

    contains

        subroutine write_ann(name,var)
            character(len=*), intent(IN) :: name
            real(prec),       intent(IN) :: var(:,:)

            call nc_write(filename,"ann_"//name,var,dim1="xc",dim2="yc",dim3="time", &
                          start=[1,1,1],count=[nx,ny,1],ncid=ncid)
        end subroutine write_ann

    end subroutine smbpal_restart_write

    subroutine smbpal_restart_read(smb,filename,ann_restored)
        ! Restore the prognostic snowpack state (H_snow, alb_s) and the annual
        ! fields of the last update. ann_restored is false for a file without
        ! the annual fields (written before they were added): smb%ann is then
        ! unset and the host must rebuild it with an update, which steps the
        ! snowpack one extra year (the restart is not exact).

        implicit none

        type(smbpal_class), intent(INOUT) :: smb
        character(len=*),   intent(IN)    :: filename
        logical,            intent(OUT)   :: ann_restored

        ! Local variables
        integer :: nx, ny

        nx = size(smb%now%H_snow,1)
        ny = size(smb%now%H_snow,2)

        call nc_read(filename,"H_snow",smb%now%H_snow,start=[1,1,1],count=[nx,ny,1])
        call nc_read(filename,"alb_s", smb%now%alb_s, start=[1,1,1],count=[nx,ny,1])

        ann_restored = nc_exists_var(filename,"ann_smb")
        if (ann_restored) then
            call read_ann("t2m",      smb%ann%t2m)
            call read_ann("pr",       smb%ann%pr)
            call read_ann("sf",       smb%ann%sf)
            call read_ann("S",        smb%ann%S)
            call read_ann("sigma",    smb%ann%sigma)
            call read_ann("PDDs",     smb%ann%PDDs)
            call read_ann("tsrf",     smb%ann%tsrf)
            call read_ann("H_snow",   smb%ann%H_snow)
            call read_ann("alb_s",    smb%ann%alb_s)
            call read_ann("smbi",     smb%ann%smbi)
            call read_ann("smb",      smb%ann%smb)
            call read_ann("melt",     smb%ann%melt)
            call read_ann("runoff",   smb%ann%runoff)
            call read_ann("refrz",    smb%ann%refrz)
            call read_ann("melt_net", smb%ann%melt_net)
        else
            write(*,*) "smbpal_restart_read:: warning: "//trim(filename)//" has no annual fields &
                       &(ann_*); they are rebuilt by an update, so the restart is not exact."
        end if

        write(*,*) "smbpal_restart_read:: read "//trim(filename)

        return

    contains

        subroutine read_ann(name,var)
            character(len=*), intent(IN)    :: name
            real(prec),       intent(INOUT) :: var(:,:)

            call nc_read(filename,"ann_"//name,var,start=[1,1,1],count=[nx,ny,1])
        end subroutine read_ann

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


