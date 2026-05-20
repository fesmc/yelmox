module pico
    ! PICO module

    use nml 
    use ncio 
    use pico_geometry
    use pico_physics

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp  = sp 

    ! Physical constants 
    real(wp), parameter :: rho_ice =  917.d0        ! Density ice           [kg/m^3] 
    real(wp), parameter :: rho_w   = 1000.d0        ! Density water         [kg/m^3] 
    real(wp), parameter :: rho_sw  = 1028.d0        ! Density seawater      [kg/m^3] 
    real(wp), parameter :: g       = 9.81d0         ! Gravitational accel.  [m/s^2]
    real(wp), parameter :: cp_o    = 3974.d0        ! Specific heat capacity of ocean mixed layer [J/kg*ºC]
    real(wp), parameter :: L_ice   = 3.34e5         ! Latent heat of fusion of ice [J/kg] 
 
    real(wp), parameter :: year_to_sec = 365.0*24.0*60.0*60.0
    real(wp), parameter :: lambda      = L_ice/cp_o ! K
    real(wp), parameter :: rho_ice_sw  = rho_ice / rho_sw 

    type pico_param_class

        real(wp) :: n_box
        real(wp) :: a_pico, b_pico, c_pico
        real(wp) :: alpha_pico, beta_pico
        real(wp) :: rho_star
        real(wp) :: gamma_tstar
        real(wp) :: C_over

        character(len=512) :: domain

    end type

    type pico_state_class 
        
        real(wp), allocatable   :: d_shlf(:,:)          ! Shelf distance to grl
        real(wp), allocatable   :: d_if(:,:)            ! Shelf distance to ice front
        real(wp), allocatable   :: r_shlf(:,:)          ! Relative ice shelf position
        real(wp), allocatable   :: n_shlf(:,:)          ! Max number of boxes per basin
        real(wp), allocatable   :: A_box(:,:)           ! Area of box
        real(wp), allocatable   :: boxes(:,:)           ! Box number of ice shelf
        real(wp), allocatable   :: bmb_shlf(:,:)        ! Shelf basal mass balance [m/a]
        real(wp), allocatable   :: T_box(:,:)           ! Temperature of box [ºC]
        real(wp), allocatable   :: S_box(:,:)           ! Salinity of box [PSU]
        real(wp), allocatable   :: CC(:,:)              ! Overtuning strength
       
    end type 

    type pico_class
        type(pico_param_class) :: par 
        type(pico_state_class) :: now 
    end type

    private
    public :: pico_class
    public :: pico_update
    public :: pico_init
    public :: pico_end

contains
    
    subroutine pico_update(pico,to,so,H_ice,z_bed,f_grnd,z_sl,basins,mask_ocn,dx)
        
        implicit none
        
        type(pico_class), intent(INOUT) :: pico
        real(wp), intent(IN) :: to(:,:)
        real(wp), intent(IN) :: so(:,:)
        real(wp), intent(IN) :: H_ice(:,:) 
        real(wp), intent(IN) :: z_bed(:,:) 
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: z_sl(:,:) 
        real(wp), intent(IN) :: basins(:,:)
        integer,  intent(IN) :: mask_ocn(:,:) 
        real(wp), intent(IN) :: dx

        ! Local variables
        integer :: m, i, j, k, l, nx, ny, ngr

        real(wp) :: bmb_floating, pm_point, pm_point_box0
        real(wp), allocatable :: T_star(:,:)
        logical,  allocatable :: is_grline(:,:)

        nx = size(f_grnd,1)
        ny = size(f_grnd,2) 

        allocate(is_grline(nx,ny))
        allocate(T_star(nx,ny))

        ! First, update the geometry 
        call pico_update_geometry(pico,H_ice,f_grnd,basins,dx)

        ! Determine location of grounding line 
        is_grline = calc_grline(f_grnd,H_ice)
 
        ! Initialize Tstar, Tbox, Sbox
        T_star = 0.0
        pico%now%T_box  = to
        pico%now%S_box  = so

        ! PISM condition in PICO. Box 0 temperature cannot be below the pressure-melting point.
        do l = 1, ny
        do k = 1, nx
            pm_point_box0 = calc_theta_pm(pico%now%S_box(k,l),pico%par%a_pico,pico%par%b_pico,pico%par%c_pico,H_ice(k,l))
            if (pico%now%T_box(k,l) .lt. pm_point_box0) pico%now%T_box(k,l) = pm_point_box0 + 0.001
        end do
        end do

        where(f_grnd .eq. 1.0)
            pico%now%T_box = 0.0
            pico%now%S_box = 0.0
        end where

        ! 1. Calculate current ice shelf bmb field (grounded-ice bmb is
        ! calculated in ice-sheet model separately) ========
        
        do j = 1, ny
        do i = 1, nx

            if (mask_ocn(i,j) .gt. 0) then
                do m = 1, int(pico%par%n_box)
                    ! jablasco: test box 1 -> .eq. -> .ge. la m
                    if (m .eq. 1.0 .and. pico%now%boxes(i,j) .ge. 1.0) then
                        ! 1. Compute box 1
                        T_star(i,j) = calc_Tstar(to(i,j),so(i,j),H_ice(i,j),pico%par%a_pico,pico%par%b_pico,pico%par%c_pico)
                        pico%now%T_box(i,j) = calc_shelf_Tbox_1(to(i,j),so(i,j),pico%now%A_box(i,j),T_star(i,j), &
                                                                pico%par%C_over,pico%par%rho_star,pico%par%gamma_tstar, &
                                                                pico%par%alpha_pico,pico%par%beta_pico)
                        pico%now%S_box(i,j) = calc_shelf_Sbox_1(to(i,j),so(i,j),pico%now%T_box(i,j))
                        ! Compute overtuning
                        pico%now%CC(i,j) = calc_overtuning(to(i,j),so(i,j),pico%now%T_box(i,j),pico%now%S_box(i,j), &
                                                           pico%par%C_over,pico%par%rho_star,pico%par%alpha_pico,pico%par%beta_pico)             

                    ! jablasco:test box 1. 1.0 -> n_box
                    else if (m .gt. 1.0 .and. pico%now%boxes(i,j) .ge. m) then
                        ! Compute rest boxes
                        T_star(i,j) = calc_Tstar(pico%now%T_box(i,j),pico%now%S_box(i,j),H_ice(i,j), &
                                                 pico%par%a_pico,pico%par%b_pico,pico%par%c_pico)
                        call calc_shelf_TS_box_n(pico%now%T_box(i,j), pico%now%S_box(i,j), pico%now%A_box(i,j), &
                                                 T_star(i,j), pico%now%CC(i,j), pico%par%a_pico, pico%par%gamma_tstar)

                    end if

                end do

                ! Compute melting
                pm_point = calc_theta_pm(pico%now%S_box(i,j),pico%par%a_pico,pico%par%b_pico,pico%par%c_pico,H_ice(i,j))
                bmb_floating = calc_melt_rate_pico(pico%now%T_box(i,j),pm_point,pico%par%gamma_tstar)
               
                ! Apply melting to purely floating points or ocean
                if(f_grnd(i,j) .eq. 0.0) pico%now%bmb_shlf(i,j) = bmb_floating

                ! Ensure that ice accretion only occurs where ice exists 
                if (pico%now%bmb_shlf(i,j) .gt. 0.0 .and. H_ice(i,j) .eq. 0.0) pico%now%bmb_shlf(i,j) = 0.0 

                ! No accretion allowed in  Box 1
                if (pico%now%bmb_shlf(i,j) .gt. 0.0 .and. pico%now%boxes(i,j) .eq. 1.0) pico%now%bmb_shlf(i,j) = 0.0

            else 
                ! Grounded point, or floating point not connected to the ocean 
                ! Set bmb_shlf to zero 

                pico%now%bmb_shlf(i,j) = 0.0
                
            end if 
 
        end do
        end do  

        return
        
    end subroutine pico_update

    subroutine pico_update_geometry(pico,H_ice,f_grnd,basins,dx)

        implicit none

        type(pico_class), intent(INOUT) :: pico
        real(wp), intent(IN) :: H_ice(:,:)
        real(wp), intent(IN) :: f_grnd(:,:)
        real(wp), intent(IN) :: basins(:,:)
        real(wp), intent(IN) :: dx

        ! Local variables
        logical,  allocatable :: is_grline(:,:), is_margin(:,:)
        real(wp), allocatable :: boxes_basin(:,:)
        real(wp) :: d_max, d_max_basin
        integer  :: m

        allocate(is_grline(size(f_grnd,1),size(f_grnd,2)))
        allocate(is_margin(size(f_grnd,1),size(f_grnd,2)))
        allocate(boxes_basin(size(f_grnd,1),size(f_grnd,2)))
        is_grline = calc_grline(f_grnd,H_ice)
        is_margin = calc_margin(H_ice)
        d_max = 0.0
        d_max_basin = 0.0

        ! compute ice shelf distance to margin and grounding-line
        call pico_calc_shelf_extent(pico%now%d_shlf,pico%now%d_if,pico%now%r_shlf,H_ice,f_grnd,is_grline,is_margin,dx)
        d_max = maxval(pico%now%d_shlf)

        do m = 1 , int(maxval(basins))
            d_max_basin = maxval(pico%now%d_shlf,basins .eq. m)
            boxes_basin = 0.0
            where(basins .eq. m)
                ! Compute boxes and area
                ! jablasco: create first indep. box basin mask seems to work
                boxes_basin = pico_calc_shelf_boxes(H_ice,f_grnd,pico%now%d_shlf,pico%now%r_shlf,pico%par%n_box,d_max,d_max_basin)
                pico%now%boxes = boxes_basin
                pico%now%A_box = pico_calc_area_boxes(boxes_basin,f_grnd,pico%par%n_box,dx)
            end where
            d_max_basin = 0.0
        end do        

    end subroutine pico_update_geometry

    ! =================
    ! 
    ! Accounting ...
    ! 
    ! =================

    subroutine pico_init(pico,filename,nx,ny,domain)
        ! Initialize pico object 

        implicit none

        type(pico_class), intent(OUT)     :: pico
        character(len=*), intent(IN)      :: filename
        integer,          intent(IN)      :: nx, ny 
        character(len=*), intent(IN)      :: domain

        ! Load parameters
        call pico_par_load(pico%par,domain,filename)

        ! Allocate the object 
        call pico_allocate(pico%now,nx,ny)

        ! Initialize variables 
        pico%now%bmb_shlf      = 0.0  
        pico%now%T_box         = 0.0
        pico%now%S_box         = 0.0

        return

    end subroutine pico_init
    
    subroutine pico_par_load(par,domain,filename,init)
        ! Load PICO parameters

        type(pico_param_class), intent(OUT) :: par
        character(len=*), intent(IN) :: domain,filename
        logical, optional :: init
        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE.

        call nml_read(filename,"pico","n_box",       par%n_box,       init=init_pars)
        call nml_read(filename,"pico","a_pico",      par%a_pico,      init=init_pars)
        call nml_read(filename,"pico","b_pico",      par%b_pico,      init=init_pars)
        call nml_read(filename,"pico","c_pico",      par%c_pico,      init=init_pars)
        call nml_read(filename,"pico","alpha_pico",  par%alpha_pico,  init=init_pars)
        call nml_read(filename,"pico","beta_pico",   par%beta_pico,   init=init_pars)
        call nml_read(filename,"pico","rho_star",    par%rho_star,    init=init_pars)
        call nml_read(filename,"pico","gamma_tstar", par%gamma_tstar, init=init_pars)
        call nml_read(filename,"pico","C_over",      par%C_over,      init=init_pars)

        ! Determine derived parameters
        par%domain = trim(domain)

        ! convert T_fp from C to K
        par%b_pico     = par%b_pico + 273.15

        ! 1Sv = 10^6 m3/s. Then 1/s to 1/yr
        par%C_over     = par%C_over*year_to_sec*10**6

        return

    end subroutine pico_par_load

    subroutine pico_end(pico)

        implicit none 

        type(pico_class) :: pico 

        ! Deallocate pico state object
        call pico_deallocate(pico%now)

        return 

    end subroutine pico_end

    subroutine pico_allocate(now,nx,ny)

        implicit none 

        type(pico_state_class) :: now 
        integer :: nx, ny 

        ! Make object is deallocated
        call pico_deallocate(now)

        ! Allocate pico 
        allocate(now%T_box(nx,ny))
        allocate(now%S_box(nx,ny))
        allocate(now%d_shlf(nx,ny))
        allocate(now%n_shlf(nx,ny))
        allocate(now%d_if(nx,ny))
        allocate(now%r_shlf(nx,ny))
        allocate(now%A_box(nx,ny))
        allocate(now%bmb_shlf(nx,ny))
        allocate(now%CC(nx,ny))
        allocate(now%boxes(nx,ny))

        ! Initialize variables 
        now%T_box    = 0.0
        now%S_box    = 0.0
        now%d_shlf   = 0.0
        now%d_if     = 0.0
        now%r_shlf   = 0.0
        now%n_shlf   = 0.0
        now%A_box    = 0.0
        now%bmb_shlf = 0.0
        now%CC       = 0.0
        now%boxes    = 0.0

        return

    end subroutine pico_allocate

    subroutine pico_deallocate(now)

        implicit none 

        type(pico_state_class) :: now 

        ! Allocate state objects
        if (allocated(now%T_box))    deallocate(now%T_box)
        if (allocated(now%S_box))    deallocate(now%S_box)
        if (allocated(now%d_shlf))   deallocate(now%d_shlf)
        if (allocated(now%d_if))     deallocate(now%d_if)
        if (allocated(now%r_shlf))   deallocate(now%r_shlf)
        if (allocated(now%n_shlf))   deallocate(now%n_shlf)
        if (allocated(now%A_box))    deallocate(now%A_box)
        if (allocated(now%bmb_shlf)) deallocate(now%bmb_shlf)
        if (allocated(now%CC))       deallocate(now%CC)
        if (allocated(now%boxes))    deallocate(now%boxes)        

        return

    end subroutine pico_deallocate
    
end module pico
