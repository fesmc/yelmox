

module basal_hydrology 
    ! This module manages the variables associated with the calculation
    ! of basal hydrology (water pressure) in the model 

    use nml 
    !use yelmo_defs !, only :: sp, dp, prec 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp 

    real(prec), parameter :: sec_year = 365.0*24.0*60.0*60.0   ! [s/a]
    real(prec), parameter :: G         = 9.81                  ! [m/s^2]
    real(prec), parameter :: pi        = 3.14159265359

    real(prec), parameter :: rho_ice   = 910.0                 ! [kg/m^3]
    real(prec), parameter :: rho_ocean = 1028.0                ! [kg/m^3] 
    real(prec), parameter :: rho_water = 1000.0                ! [kg/m^3] 

    real(prec), parameter :: rho_ice_G   = rho_ice*G           ! [kg/m^3 (m/s^2)]
    real(prec), parameter :: rho_ocean_G = rho_ocean*G         ! [kg/m^3 (m/s^2)]
    real(prec), parameter :: rho_water_G = rho_water*G         ! [kg/m^3 (m/s^2)]

    
    type hydro_param_class

        integer :: init_method
        real(prec) :: H0
        real(prec) :: H0_ice_min, H0_ice_max 
        real(prec) :: H0_min, H0_max 

        integer :: method
        real(prec) :: H_water_max      ! Maximum allowed water thickness [m]      
        real(prec) :: till_max         ! Maximum till thickness [m]
        real(prec) :: till_porosity
        real(prec) :: till_infiltr
        real(prec) :: kond            ! Initial conductivity [m/s]

        ! Internal parameters 
        real(prec) :: till_Hwat_max    ! Maximum water allowed in sediment       
        real(prec) :: kond_max
        real(prec) :: keff_max   

        real(prec) :: H_buffer_max     ! Maximum water thickness in the buffer [m]
    end type 

    type hydro_state_class
        real(prec) :: time, dt  
        real(prec), allocatable :: H_water(:,:)        ! Water pressure [m water equivalent]
        real(prec), allocatable :: Hdot_water(:,:)     ! Rate of change in water pressure [m/a]
        real(prec), allocatable :: p_water(:,:)        ! Water pressure [Pa]
        real(prec), allocatable :: hw(:,:)             ! Water depth in the till [m]

        real(prec), allocatable :: H_water_buffer(:,:) ! Buffer to avoid oscillations at low water values
    end type 

    type hydro_class
        type(hydro_param_class) :: par 
        type(hydro_state_class) :: now 
    end type

    private
    public :: hydro_class
    public :: hydro_init 
    public :: hydro_init_state 
    public :: hydro_update
!     public :: hydro_end 
    
    public :: relaxation_waterdif





contains 

    subroutine hydro_update(hyd,H_ice,z_bed,z_srf,z_sl,bmb_ice,f_grnd,dx,dy,time)
        ! Update the state of the hydro variables 

        implicit none 

        type(hydro_class), intent(INOUT) :: hyd 
        real(prec), intent(IN) :: H_ice(:,:)      ! [m] Ice thickness 
        real(prec), intent(IN) :: z_bed(:,:)      ! [m] Bedrock elevation
        real(prec), intent(IN) :: z_srf(:,:)      ! [m] Surface elevation
        real(prec), intent(IN) :: z_sl(:,:)       ! [m] Sea level 
        real(prec), intent(IN) :: bmb_ice(:,:)    ! [m/a] Basal mass balance of ice
        real(prec), intent(IN) :: f_grnd(:,:)     ! [-] Grounded fraction of grid cell
        real(prec), intent(IN) :: dx, dy          ! [m], [m] Grid cell resolutions
        real(prec), intent(IN) :: time            ! [a] Current external time to advance to

        ! Local variables
        integer :: i, j, nx, ny  
        real(prec), allocatable :: z_base(:,:)     ! [m] z_srf - H_ice 
        real(prec), allocatable :: H_ocean(:,:)    ! [m] ice-free: sealevel-z_bed, ice-covered: z_base-z_bed
        logical,    allocatable :: is_float(:,:)   ! Floating?
        logical,    allocatable :: is_grz(:,:)     ! Grounding zone?
        real(prec), allocatable :: bmb_w(:,:)      ! [m/a] Basal mass balance of water
        
        nx = size(H_ice,1)
        ny = size(H_ice,2) 

        ! Allocate variables 
        allocate(z_base(nx,ny))
        allocate(H_ocean(nx,ny))
        allocate(is_float(nx,ny))
        allocate(is_grz(nx,ny))
        allocate(bmb_w(nx,ny))

        ! Determine current time step and time 
        hyd%now%dt   = max(time - hyd%now%time,0.0)
        hyd%now%time = time 

        ! Determine basal water mass balance 
        bmb_w   = -bmb_ice*rho_water/rho_ice   
        
        ! Determine ice base 
        z_base = z_srf - H_ice 

        if (hyd%now%dt .gt. 0.0) then 

            ! Calculate floating mask 
            is_float = f_grnd .eq. 0.0 

            ! Calculate ocean depth
            where (.not. is_float)
                ! Grounded land point, no ocean
                H_ocean = 0.0 
            else where (H_ice .gt. 0.0)
                ! Floating ice point
                H_ocean = z_base - z_bed
            elsewhere 
                ! Opean ocean point 
                H_ocean = z_sl - z_bed 
            end where  

            ! Make sure H_ocean makes sense 
            where( H_ocean .lt. 0.0) H_ocean = 0.0 

            is_grz = .FALSE. 
            do i = 2, nx-1
            do j = 2, ny-1 
                if (.not. is_float(i,j) .and. &
                    (is_float(i-1,j).or.is_float(i+1,j).or. &
                     is_float(i,j-1).or.is_float(i,j+1))) then 
                    is_grz(i,j) = .TRUE. 
                end if 

            end do 
            end do 

            ! Update the hydro object depending on the method desired
            select case(hyd%par%method)

                case(0)
                    ! Constant H_water at initial value 
                    
                    ! Pass, nothing happens

                case(1)
                    ! Local mass balance of H_water 

                    call calc_basal_water_local(hyd%par,hyd%now%H_water,H_ice,H_ocean,bmb_w,is_float,hyd%now%dt)

                case(2)
                    ! 2D Diffusion of H_water 

                    call calc_basal_water_diffusion(hyd%par,hyd%now%H_water,hyd%now%Hdot_water, &
                                    hyd%now%p_water,hyd%now%hw,H_ice,H_ocean,z_base,z_srf,z_sl,bmb_w, &
                                    is_float,is_grz,dx,dy,hyd%now%dt)

                case DEFAULT 

                    write(*,*) "hydro_update:: error: method must be in (0,1,2)."
                    write(*,*) "method = ", hyd%par%method 
                    stop 

            end select 

!             ! Adjust water thickness based on buffer 
!             where(hyd%now%H_water_buffer .gt. hyd%par%H_buffer_max)

!                 hyd%now%H_water =hyd%now%H_water_buffer

!             elsewhere 
!                 hyd%now%H_water = 0.0 

!             end where 
            
!             where(hyd%now%H_water .lt. hyd%par%H_buffer_max)
!                 hyd%now%H_water = 0.0 
!             end where 
        
        end if 


        return 

    end subroutine hydro_update 

    subroutine hydro_init(hyd,filename,nx,ny)

        implicit none 

        type(hydro_class), intent(INOUT) :: hyd 
        character(len=*),  intent(IN)    :: filename 
        integer,           intent(IN)    :: nx, ny 

        ! Load parameter options 
        call hydro_par_load(hyd%par,filename)

        ! Initialize state object 
        call hydro_allocate(hyd%now,nx,ny)

        ! Set fields to zero for now
        ! (use hydro_init_state later to intialize properly)
        hyd%now%H_water        = 0.0 
        hyd%now%Hdot_water     = 0.0 
        hyd%now%p_water        = 0.0 
        hyd%now%hw             = 0.0 
        
        hyd%now%H_water_buffer = 0.0

        hyd%now%time = 1e10     ! Set time way into the future, for now 
        hyd%now%dt   = 0.0 

        return 

    end subroutine hydro_init 

    subroutine hydro_init_state(hyd,H0,H_ice,H_ocean,time)
        ! Initialize the state of the hydro fields
        ! (only after calling hydro_init)

        implicit none 

        type(hydro_class), intent(INOUT) :: hyd 
        real(prec), intent(IN), optional :: H0(:,:)
        real(prec), intent(IN), optional :: H_ice(:,:), H_ocean(:,:)  
        real(prec), intent(IN) :: time 

        if (present(H0)) then 
            ! Initialize using the provided field

            ! Consistency check 
            if (size(H0,1) .ne. size(hyd%now%H_water,1) .or.  & 
                size(H0,2) .ne. size(hyd%now%H_water,2)) then 

                write(*,*) "hydro_init_state:: error: H0 field dimensions must match the &
                           &hydro array dimensions."
                write(*,*) "dim(H0):    ", size(H0,1), size(H0,2)
                write(*,*) "dim(hydro): ", size(hyd%now%H_water,1), size(hyd%now%H_water,2)
                stop 
            end if 
            
            ! Set H_water equal to H0
            hyd%now%H_water = H0 

        else 
            ! Initialize using the method of choice (par%init_method)

            select case(hyd%par%init_method)

                case(0)
                    ! Intially zero everywhere
                    hyd%now%H_water = 0.0 

                case(1)
                    ! Spatially constant intial value 
                    hyd%now%H_water = hyd%par%H0

                case(2)
                    ! H_water = f(H_ice)

                    ! Consistency checks
                    if (.not. present(H_ice)) then 
                        write(*,*) "hydro_init_state:: error: for init_method=2, &
                                   &the ice thickness field argument 'H_ice' must be present." 
                        stop 
                    end if 

                    if (.not. present(H_ocean)) then 
                        write(*,*) "hydro_init_state:: error: for init_method=2, &
                                   &the ocean depth field argument 'H_ocean' must be present." 
                        stop 
                    end if 

                    if (size(H_ice,1) .ne. size(hyd%now%H_water,1) .or.  & 
                        size(H_ice,2) .ne. size(hyd%now%H_water,2)) then 

                        write(*,*) "hydro_init_state:: error: H_ice field dimensions must match the &
                                   &hydro array dimensions."
                        write(*,*) "dim(H_ice): ", size(H_ice,1), size(H_ice,2)
                        write(*,*) "dim(hydro): ", size(hyd%now%H_water,1), size(hyd%now%H_water,2)
                        stop 
                    end if 
                
                    ! H_water is zero everywhere, except it is linearly
                    ! dencreasing from H0_max to H0_min for ice thicknesses in 
                    ! range of H0_ice_min to H0_ice_max. 
                    hyd%now%H_water = 0.0 

                    where(H_ice .lt. hyd%par%H0_ice_max .and. &
                          H_ice .gt. hyd%par%H0_ice_min)

                        hyd%now%H_water = hyd%par%H0_max + (hyd%par%H0_min-hyd%par%H0_max)  &
                                            *(H_ice-hyd%par%H0_ice_min) / (hyd%par%H0_ice_max-hyd%par%H0_ice_min)

                    end where

                case DEFAULT 

                    write(*,*) "hydro_init_state:: error: initialization method must be one of (0,1,2)."
                    write(*,*) "init_method = ", hyd%par%init_method 
                    stop

            end select

            ! Handle floating cases 
            where (H_ocean*rho_ocean/rho_water - H_ice .gt. 0.0) 
                hyd%now%H_water = H_ocean 
            end where 

        end if 

        ! Set the current time and time step 
        hyd%now%time = time 
        hyd%now%dt   = 0.0 

        return 

    end subroutine hydro_init_state 


    subroutine hydro_par_load(par,filename,init)

        implicit none 

        type(hydro_param_class)  :: par
        character(len=*)         :: filename 
        logical, optional        :: init

        ! Local variables 
        logical :: init_pars 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE. 
        
        call nml_read(filename,"basal_hydrology","init_method",   par%init_method,  init=init_pars)
        call nml_read(filename,"basal_hydrology","H0",            par%H0,           init=init_pars)
        call nml_read(filename,"basal_hydrology","H0_ice_min",    par%H0_ice_min,   init=init_pars)
        call nml_read(filename,"basal_hydrology","H0_ice_max",    par%H0_ice_max,   init=init_pars)
        call nml_read(filename,"basal_hydrology","H0_min",        par%H0_min,       init=init_pars)
        call nml_read(filename,"basal_hydrology","H0_max",        par%H0_max,       init=init_pars)
        call nml_read(filename,"basal_hydrology","method",        par%method,       init=init_pars)
        call nml_read(filename,"basal_hydrology","H_water_max",   par%H_water_max,  init=init_pars)
        call nml_read(filename,"basal_hydrology","till_max",      par%till_max,     init=init_pars)
        call nml_read(filename,"basal_hydrology","till_porosity", par%till_porosity,init=init_pars)
        call nml_read(filename,"basal_hydrology","till_infiltr",  par%till_infiltr, init=init_pars)
        call nml_read(filename,"basal_hydrology","kond",          par%kond,         init=init_pars)

        ! Set some additional internal parameters 

        ! Maximum height the water can reach in the till layer
        par%till_Hwat_max = par%till_max*par%till_porosity   ! [m]

        ! Maximum conductivities
        par%kond_max    = 1.0*sec_year                    ! [m/a]
        par%keff_max    = par%kond_max*par%till_Hwat_max  ! [m^2/a]

        par%H_buffer_max = 0.1  ! [m]

        return 

    end subroutine hydro_par_load

    subroutine hydro_allocate(now,nx,ny)

        implicit none 

        type(hydro_state_class), intent(INOUT) :: now 
        integer :: nx, ny 

        ! First make sure fields are deallocated
        call hydro_deallocate(now)

        ! Allocate fields to desired dimensions
        allocate(now%H_water(nx,ny))
        allocate(now%Hdot_water(nx,ny))
        allocate(now%p_water(nx,ny))
        allocate(now%hw(nx,ny))

        allocate(now%H_water_buffer(nx,ny))
        
        ! Initialize all fields to zero
        now%H_water        = 0.0
        now%Hdot_water     = 0.0
        now%p_water        = 0.0
        now%hw             = 0.0
        now%H_water_buffer = 0.0
        
        return 

    end subroutine hydro_allocate 

    subroutine hydro_deallocate(now)

        implicit none 

        type(hydro_state_class), intent(INOUT) :: now 

        if (allocated(now%H_water))    deallocate(now%H_water)
        if (allocated(now%Hdot_water)) deallocate(now%Hdot_water)
        if (allocated(now%p_water))    deallocate(now%p_water)
        if (allocated(now%hw))         deallocate(now%hw)

        if (allocated(now%H_water_buffer)) deallocate(now%H_water_buffer)
        
        return 

    end subroutine hydro_deallocate 



    ! ==== BASAL HYDROLOGY PHYSICS ===============================

    subroutine calc_basal_water_local(par,H_water,H_ice,H_ocean,bmb_w,is_float,dt)
        ! method=1: local mass balance of H_water
        ! Formerly inside of `eaubasale`
        ! Note: only call when ISYNCHRO==1

        implicit none 

        type(hydro_param_class), intent(IN) :: par 
        real(prec), intent(INOUT) :: H_water(:,:)
        real(prec), intent(IN)    :: H_ice(:,:)
        real(prec), intent(IN)    :: H_ocean(:,:) 
        real(prec), intent(IN)    :: bmb_w(:,:)       ! [m/a] basal water mass balance
        logical, intent(IN)    :: is_float(:,:) 
        real(prec), intent(IN)    :: dt 

        where (.not. is_float .and. H_ice .ge. 0.1)
            ! Grounded ice point

            ! Update mass balance of H_water
            H_water = H_water + dt*bmb_w - dt*par%till_infiltr

            ! Restrict H_water to values within limits
            H_water = max(H_water,0.0)
            H_water = min(H_water,par%H_water_max)

        else where (.not. is_float) 
            ! Ice-free land above sea level 

            H_water = 0.0 

        elsewhere
            ! Set water pressure to maximum (ie, ocean depth)

            H_water = H_ocean

        end where 

        return 

    end subroutine calc_basal_water_local

    subroutine calc_basal_water_diffusion(par,H_water,Hdot_water,p_water,hw, &
                                H_ice,H_ocean,z_base,z_srf,z_sl,bmb_w,is_float,is_grz,dx,dy,dt)
        ! version correspondant à la thèse de Vincent
        ! (previously `eaubasale` in eaubasale-0.5_mod.f90)
        ! Note: only call when nt>0 and dt>0

        implicit none 

        type(hydro_param_class), intent(IN) :: par 
        real(prec), intent(INOUT) :: H_water(:,:)    ! current water pressure
        real(prec), intent(INOUT) :: Hdot_water(:,:) ! Rate of change [m/a]
        real(prec), intent(INOUT) :: p_water(:,:)    ! water pressure [Pa]
        real(prec), intent(INOUT) :: hw(:,:)         ! till infiltration 
        
        real(prec), intent(IN) :: H_ice(:,:)      ! ice thickness
        real(prec), intent(IN) :: H_ocean(:,:)    ! sealevel-Bsoc
        real(prec), intent(IN) :: z_base(:,:)      ! B
        real(prec), intent(IN) :: z_srf(:,:)      ! S [m]
        real(prec), intent(IN) :: z_sl(:,:)            ! sealevel [m]
        real(prec), intent(IN) :: bmb_w(:,:)        ! [m/a] basal water mass balance
!         real(prec), intent(IN) :: uvb(:,:)        ! Basal velocity magnitude 
        logical, intent(IN) :: is_float(:,:)     ! Floating
        logical, intent(IN) :: is_grz(:,:)       ! Grounding zone 
        real(prec), intent(IN) :: dx, dy          ! Grid resolution [m] 
        real(prec), intent(IN) :: dt              ! time step 

        ! Local variables
        integer :: nx, ny, i, j 
        real(prec) :: neff_tmp, z_float  

        integer, allocatable :: klimit(:,:) 
        real(prec), allocatable :: limit_hw(:,:)  
        real(prec), allocatable :: kond(:,:)
        real(prec), allocatable :: keff(:,:)
        real(prec), allocatable :: phiwx(:,:), phiwy(:,:)
        real(prec), allocatable :: pgx(:,:), pgy(:,:) 
        real(prec), allocatable :: H_water_old(:,:) 
        real(prec), allocatable :: pot_f(:,:) 
        real(prec), allocatable :: pot_w(:,:) 
        logical, allocatable :: is_float_mx(:,:), is_float_my(:,:)  
        
        ! Allocate local variables 
        nx = size(H_water,1)
        ny = size(H_water,2)

        allocate(klimit(nx,ny),limit_hw(nx,ny),kond(nx,ny),keff(nx,ny), &
                 phiwx(nx,ny),phiwy(nx,ny),pgx(nx,ny),pgy(nx,ny),H_water_old(nx,ny), &
                 pot_f(nx,ny),pot_w(nx,ny))
        allocate(is_float_mx(nx,ny),is_float_my(nx,ny)) 

        ! Set initial values 
        H_water_old = H_water
        kond        = par%kond*sec_year

!         ! Adjust conductivity mask as a function of basal velocity
!         if par%

        ! For now, treat x/y mask the same 
        is_float_mx = is_float 
        is_float_my = is_float 

        ! Initialize all other new fields to zero 
        klimit   = 0
        limit_hw = 0.0  
        hw       = 0.0 
        keff     = 0.0 
        phiwx    = 0.0 
        phiwy    = 0.0 
        pgx      = 0.0 
        pgy      = 0.0 
        pot_f    = 0.0 
        pot_w    = 0.0 

        where (is_float)  
            ! Floating points
            klimit   = 1
            limit_hw = H_ocean*(rho_ocean_G)/(rho_water_G)

        else where((.not.is_float).and.(H_ice.lt.1.0))  
            ! Ice margin
            klimit   = 1
            limit_hw = 10.0    ! River of depth=10m
        
        elsewhere 

            ! No ice 
            klimit   = 0 
            limit_hw = -9999.0 

        end where 

        ! Calculate boundary conditions
        H_water = max(H_water,0.0)
        hw      = min(H_water,par%till_Hwat_max)

        ! Determine boundary pressure 

        ! Gravitational potential plus pressure from ice (to force diffusion)
        pot_w = rho_water_G*z_base + rho_ice_G*H_ice
     
        ! Pressure at the base of the ice shelf   
        pot_f = rho_ocean_G*(z_sl-z_srf+H_ice)  ! pression a la base de l'ice shelf
!         pot_f = rho_ocean_G* (0.910*H_ice)   ! ajr: to avoid passing z_sl and z_srf, should be the same right??
 
        ! Grisli global equivalents
        ! z_sl   = sealevel
        ! z_srf  = S 
        ! z_base = B
        ! H_ice  = H 
        
        ! Calculate the pressure gradient 
        do j = 2, ny
        do i = 2, nx

            if (H_ice(i,j).gt.25.0) then
                ! ajr: hard-coded threshold - only calculate if H_ice > 25 m... why?? 

                ! x-direction
                if (is_float_mx(i,j)) then 
                    pgx(i,j)=(pot_f(i-1,j)-pot_f(i,j))/dx
                else
                    pgx(i,j)=(pot_w(i-1,j)-pot_w(i,j))/dx
                endif

                ! y-direction 
                if (is_float_my(i,j)) then 
                    pgy(i,j)=(pot_f(i,j-1)-pot_f(i,j))/dy
                else
                    pgy(i,j)=(pot_w(i,j-1)-pot_w(i,j))/dy
                endif

            else 
                pgx(i,j) = 0.0 
                pgy(i,j) = 0.0 

            endif

            ! Make pgx/pgy unitless 
            pgx(i,j) = pgx(i,j) / rho_water_G
            pgy(i,j) = pgy(i,j) / rho_water_G

        end do
        end do

        if (dt.gt.0.0) then  !!!!!!!!!!!!!!!!!!!!!!relax_water si dt>0

        ! 1. Determine conductivity of point =========================
        ! Note: grounding line points have high hydraulic conductivity even if the base is cold
        ! (note from eaubasale-0.5_mod.F90: modif catritz 18 janvier 2005)

        do j = 2, ny-1
        do i = 2, nx-1

            ! Maximum (ie, infinite) conductivity with small ice thickness or floating point
            if (is_float(i,j).or.H_ice(i,j).le.1.5)  kond(i,j) = par%kond_max
        
            ! Also high conductivity in the grounding zone 
            if (is_grz(i,j)) kond(i,j) = par%kond_max 

            ! High conductivity when effective pressure is small (~100 bar)
            neff_tmp = 0.91*H_ice(i,j)-H_water(i,j)
            neff_tmp = max(100.0,neff_tmp)
            if (neff_tmp.le.1000.0) then
                kond(i,j) = kond(i,j)*1000.0/neff_tmp
            endif
            kond(i,j) = min(kond(i,j),par%kond_max)  

            ! Calculate effective conductivity (kond * pipe size, in [m^2/a])
            keff(i,j) = kond(i,j)*hw(i,j)

            ! Also high conductivity in the grounding zone 
            if (is_grz(i,j)) keff(i,j) = par%keff_max 

        end do
        end do

        ! Grid border conditions 
        kond(1,:)  = par%kond_max
        kond(nx,:) = par%kond_max
        kond(:,1)  = par%kond_max
        kond(:,ny) = par%kond_max

        ! 2. Diffuse H_water =======================================

        ! Store current H_water in old field
        H_water_old = H_water

        ! Perform diffusion
        call relaxation_waterdif(nx,ny,dt,dx,H_water_old,limit_hw,klimit,bmb_w, &
                                 par%till_infiltr,pgx,pgy,keff,par%keff_max,H_water)


        end if 
 
        ! 3. Calculate the relaxation, boundary conditions and extreme values ==========
        ! Note: Rate is calculated in 2 timesteps

        ! Initially estimate rate to be able to precisely calculate it later
        if (dt .gt. 0) then 
            Hdot_water = (H_water-H_water_old)/dt 
        end if 
                                                                                  
        ! Calculate H_water and p_water 
        do i=1,nx
        do j=1,ny

            if (is_float(i,j) .or. H_ice(i,j) .le. 1.5) then 

                if (is_float(i,j)) then
                    ! Floating ice shelves (pressure equal to water depth)
     
                    H_water(i,j) = max(H_ocean(i,j),0.0)
                    p_water(i,j) = H_water(i,j)*rho_ocean_G 
                
                else
                    ! Bare land (no water pressure)

                    H_water(i,j) = 0.0 
                    p_water(i,j) = 0.0 
                    
                end if 
                
            else
                ! Under the dome (grounded ice)

                z_float = H_ice(i,j)*rho_ice/rho_water

                if (H_water(i,j).le.0.0) then 
                    
                    H_water(i,j) = 0.0
                    p_water(i,j) = 0.0 

                else if (H_water(i,j).gt.z_float) then

                    H_water(i,j) = z_float
                    hw(i,j)      = min(H_water(i,j),par%till_Hwat_max)
                    p_water(i,j) = H_water(i,j)*rho_water_G

                end if

            end if        

            ! bloc qui pourrait servir pour mettre l'eau encore plus sous pression
            ! -----------------------------------------------------------------------------
            !  if (HWATER(i,j).gt.poro_till*hmax_till) then
            !    pwater(i,j)=pwater(i,j)+(HWATER(i,j)-poro_till*hmax_till)/(compress_w*hmax_till)
            !  endif

        end do
        end do

        ! Specify values of H_water at the corners
        H_water(1,1)   = (H_water(1,2)+H_water(2,1))/2.0
        H_water(1,ny)  = (H_water(1,ny-1)+H_water(2,ny))/2.0
        H_water(nx,1)  = (H_water(nx,2)+H_water(nx-1,1))/2.0
        H_water(nx,ny) = (H_water(nx,ny-1)+H_water(nx-1,ny))/2.0
    
        ! For water outflow
        do j=2,ny
        do i=2,nx

            ! x-direction
            if  ( abs(keff(i,j)+keff(i-1,j)) .lt. 1e-10 ) then
                phiwx(i,j) = 0.0 ! to avoid division by 0             
            else
                phiwx(i,j) = (pgx(i,j)+(H_water(i-1,j)-H_water(i,j))/dx)
                phiwx(i,j) = phiwx(i,j)*2.0*(keff(i,j)*keff(i-1,j))/(keff(i,j)+keff(i-1,j))
            endif
            pgx(i,j) = (pgx(i,j)+(H_water(i-1,j)-H_water(i,j))/dx)

            ! y-direction
            if  ( abs(keff(i,j)+keff(i,j-1)) .lt. 1e-10 ) then
                phiwy(i,j) = 0.0 ! to avoid division by 0
            else
                phiwy(i,j) = (pgy(i,j)+(H_water(i,j-1)-H_water(i,j))/dy)
                phiwy(i,j) = phiwy(i,j)*2.0*(keff(i,j)*keff(i,j-1))/(keff(i,j)+keff(i,j-1))
            endif
            pgy(i,j)=(pgy(i,j)+(H_water(i,j-1)-H_water(i,j))/dy)

        enddo
        enddo

!         ! Check boundary variables 
!         write(*,*) "calc_basal_water_diffusion:: "
!         write(*,*) "range(bmb_w):   ", minval(bmb_w), maxval(bmb_w)
!         write(*,*) "range(H_water): ", minval(H_water), maxval(H_water)
!         write(*,*) "range(p_water): ", minval(p_water), maxval(p_water)
!         write(*,*) "range(hw):      ", minval(hw), maxval(hw)
!         write(*,*) "range(tetar):   ", minval(tetar), maxval(tetar)
!         write(*,*) "range(pgx):     ", minval(pgx), maxval(pgx)
!         write(*,*) "range(pgy):     ", minval(pgy), maxval(pgy)
!         write(*,*) "range(pot_w):   ", minval(pot_w), maxval(pot_w)
!         write(*,*) "range(pot_f):   ", minval(pot_f), maxval(pot_f)
!         write(*,*) "range(pgx):     ", minval(pgx), maxval(pgx)
!         write(*,*) "range(pgy):     ", minval(pgy), maxval(pgy)
!         write(*,*) "range(kond):    ", minval(kond), maxval(kond)
!         write(*,*) "range(keff):    ", minval(keff), maxval(keff)
!         write(*,*) "range(phiwx):   ", minval(phiwx), maxval(phiwx)
!         write(*,*) "range(phiwy):   ", minval(phiwy), maxval(phiwy)
        
        return

    end subroutine calc_basal_water_diffusion

    subroutine relaxation_waterdif(NXX,NYY,DT,DX,vieuxHWATER,limit_hw,klimit, &
                                    BMB,INFILTR,PGMX,PGMY,KOND,KONDMAX,HWATER)
        ! Subroutine to relax the water pressure hwater according 
        ! to the following relationship:
        ! dhwat/dt = bmb-infiltr-d/dx(Kond*dhwat/dx)+d/dx(Kond*pgx)

        implicit none


        ! declaration des variables en entree
        !------------------------------------------------
        INTEGER, intent(in) :: NXX, NYY       ! defini la taille des tableaux
        REAL,    intent(in) ::  DT            ! pas de temps court
        REAL,    intent(in) ::  DX            ! pas en x
        REAL,    intent(in) :: INFILTR        ! basal infiltration (lose of water) 
        REAL,    intent(in) :: KONDMAX        ! maximum hydaulic conductivity (outside ice sheet)

        REAL,dimension(NXX,NYY), intent(in) :: limit_hw    ! conditions aux limites
        integer,dimension(NXX,NYY), intent(in) :: klimit    ! ou appliquer les conditions 
        REAL,dimension(NXX,NYY), intent(in) :: vieuxHWATER    ! H au pas de temps precedent 'o'
        REAL,dimension(NXX,NYY), intent(in) :: BMB      ! basal water production  'o'
        REAL,dimension(NXX,NYY), intent(in) :: PGMX     ! hydaulic potential gratient '>' 
        REAL,dimension(NXX,NYY), intent(in) :: PGMY     ! hydaulic potential gratient '^'
        REAL,dimension(NXX,NYY), intent(in) :: KOND     ! hydaulic conductivity 'o'

        ! declaration des variables en sortie
        !-------------------------------------
        REAL,dimension(NXX,NYY), intent(out):: HWATER      ! basal water thickness  'o' (pressure equivalent)


        ! declaration des variables locales
        !----------------------------------
        INTEGER :: I,J
        REAL  :: TESTH 
        REAL,dimension(NXX,NYY) :: ARELAX,BRELAX,CRELAX,DRELAX,ERELAX,FRELAX
        REAL,dimension(NXX,NYY) :: DELTAH
        REAL :: RESTE,DELH,VH
        INTEGER :: ntour
        INTEGER :: mbord
        REAL  :: DTSRGDX,dtwdx2
        LOGICAL :: STOPP
        REAL,dimension(NXX,NYY) :: KMX, KMY

        
        ! First, save previous timestep field into current one
        HWATER(:,:)= vieuxHWATER(:,:)

        ! calcul de kmx et kmx a partir de KOND
        ! conductivite hyrdraulique sur les noeuds mineurs
        ! moyenne harmonique
        ! ----------------------------------------

        do j=2,nyy
        do i=2,nxx

            if  ((kond(i,j).lt.1.e-20).or.(kond(i-1,j).lt.1.e-20)) then
                kmx(i,j)=0.0 ! to avoid division by 0
            else
                kmx(i,j)=2*(kond(i,j)*kond(i-1,j))/(kond(i,j)+kond(i-1,j))
            endif

        end do
        end do

        do j=2,nyy
        do i=2,nxx

            if  ((kond(i,j).lt.1.e-20).or.(kond(i,j-1).lt.1.e-20)) then
                kmy(i,j)=0.0 ! to avoid division by 0
            else
                kmy(i,j)=2*(kond(i,j)*kond(i,j-1))/(kond(i,j)+kond(i,j-1))
            endif

        enddo
        enddo

        ! attribution des coefficients  arelax ....
        ! ----------------------------------------
        !  SECYEAR=365.*24.*3600.
        ! rho=910.
        !  rhow=1000.
        !  rhog=rhow*9.81
        !  dtsrgdx=dt/(rhog*DX) a mon avis c'est rhow qu'il fallait utiliser. Maintenant cette 
        !  division est faite dans eaubasale 

        dtsrgdx = dt/DX
        dtwdx2  = dt/dx/dx

        arelax = 0.0
        brelax = 0.0
        crelax = 1.0
        drelax = 0.0
        erelax = 0.0
        frelax = limit_hw 
 
        do J=2,NYY-1
        do I=2,NXX-1

            if (klimit(i,j).eq.0) then
                ! Only perform caculations where water limit has not been exceeded

                ! calcul du vecteur

                FRELAX(i,j)= VIEUXHWATER(i,j)+DT*(BMB(i,j)-INFILTR)
                frelax(i,j)=frelax(i,j)+(kmx(i,j)*pgmx(i,j)-kmx(i+1,j)*pgmx(i+1,j))*dtsrgdx
                frelax(i,j)=frelax(i,j)+(kmy(i,j)*pgmy(i,j)-kmy(i,j+1)*pgmy(i,j+1))*dtsrgdx

                ! calcul des diagonales      
                arelax(i,j)=-kmx(i,j)*dtwdx2        ! arelax : diagonale i-1,j 

                brelax(i,j)=-kmx(i+1,j)*dtwdx2      ! brelax : diagonale i+1,j 

                drelax(i,j)=-kmy(i,j)*dtwdx2        ! drelax : diagonale i,j-1 

                erelax(i,j)=-kmy(i,j+1)*dtwdx2      ! drelax : diagonale i,j+1 

                crelax(i,j)=1.+((kmx(i,j)+kmx(i+1,j))+(kmy(i,j+1)+kmy(i,j+1)))*dtwdx2 
                                                    !crelax : diagonale i,j 

            else if (klimit(i,j).eq.1) then
                ! Set hwater equal to the limit in these regions

                hwater(i,j)=limit_hw(i,j)

            endif

        end do
        end do




        ! Relaxation loop :
        ! ----------------------

        ntour = 0
        stopp = .false.

        Do  while(.not.stopp)
            
            ntour=ntour+1

            do j=2,NYY-1
            do i=2,NXX-1

                RESTE = (((ARELAX(i,j)*HWATER(i-1,j) + BRELAX(i,j)*HWATER(i+1,j)) &
                      + (DRELAX(i,j)*HWATER(i,j-1) + ERELAX(i,j)*HWATER(i,j+1))) &
                      + CRELAX(i,j)*HWATER(i,j))- FRELAX(i,j)

                DELTAH(i,j) = RESTE/CRELAX(i,j)             

            end do
            end do

      

            ! il faut faire le calcul suivant dans une autre boucle car RESTE est fonction
            ! de hwater sur les points voisins.
            do j=2,NYY-1
            do i=2,NXX-1

                HWATER(i,j) = HWATER(i,j) - DELTAH(i,j)

            end do
            end do


            ! Stopping criterion:

            Delh = 0
            Vh   = 0

            DO j=2,NYY-1
            DO i=2,NXX-1

                Delh=Delh+deltah(i,j)**2

            END DO
            END DO

            if (delh.gt.0.0) then
                testh=sqrt(Delh)/((NXX-2)*(NYY-2))
            else
                testh = 0.0
            endif
            
            stopp = (testh.lt.1.e-3).or.(ntour.gt.1000)
  
        end do

        return

    end subroutine relaxation_waterdif

end module basal_hydrology 



