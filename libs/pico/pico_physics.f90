module pico_physics
    ! Module to simulate the physics via PICO

    use nml 
    use ncio 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp 

    ! Physical constants 
    real(prec), parameter :: rho_ice =  917.d0       ! Density ice           [kg/m^3] 
    real(prec), parameter :: rho_w   = 1000.d0       ! Density water         [kg/m^3] 
    real(prec), parameter :: rho_sw  = 1028.d0       ! Density seawater      [kg/m^3] 
    real(prec), parameter :: g       = 9.81d0        ! Gravitational accel.  [m/s^2]
    real(prec), parameter :: cp_o    = 3974.d0       ! Specific heat capacity of ocean mixed layer [J/kg*ºC]
    real(prec), parameter :: L_ice   = 3.34e5        ! Latent heat of fusion of ice [J/kg] 
 
    real(prec), parameter :: year_to_sec = 365.0*24.0*60.0*60.0
    real(prec), parameter :: lambda      = L_ice/cp_o ! [ºC -> K]?
    real(prec), parameter :: rho_ice_sw  = rho_ice / rho_sw  

    private
    public :: calc_Tstar
    public :: calc_shelf_Tbox_1
    public :: calc_shelf_Sbox_1
    public :: calc_shelf_TS_box_n
    public :: calc_overtuning
    public :: calc_theta_pm
    public :: calc_T_pm
    public :: apply_c_deep
    public :: calc_melt_rate_pico

contains 

    ! =======================================
    ! 
    ! Subroutines to compute pico physics 
    ! 
    ! =======================================

    function calc_Tstar(to,so,H_ice,a_pico,b_pico,c_pico) result(T_star)

        implicit none

        real(prec), intent(IN)  :: to, so, H_ice
        real(prec), intent(IN)  :: a_pico, b_pico, c_pico
        real(prec) :: T_star

        T_star = a_pico*so + b_pico - c_pico*rho_ice*g*H_ice - to
        ! Ensures that temperature input for grounding line box should not be below pressure melting point
        ! This ensures that later equations are well solvable.
        if(T_star .gt. 0.0) T_star = 0.0 !-0.0001

        return

    end function calc_Tstar

    ! First T_box1
    function calc_shelf_Tbox_1(to,so,A_box,T_star,C_over,rho_star,gamma_tstar,alpha_pico,beta_pico) result(T_box)

        implicit none

        real(prec), intent(IN)  :: A_box, T_star, to, so
        real(prec), intent(IN)  :: C_over, rho_star, gamma_tstar, alpha_pico, beta_pico
        real(prec) :: T_box

        ! Internal variables
        real(prec) :: g1,s1,p,q,D

        g1 = gamma_tstar*A_box
        s1 = so/(lambda*rho_ice_sw)
        p = g1 / (C_over * rho_star * (beta_pico * s1 - alpha_pico))
        q = p * T_star

        D = 0.25*p*p - q
        if (D .lt. 0.0) D = 0.0
      
        T_box = to - (-0.5*p + SQRT(D))

        return

    end function calc_shelf_Tbox_1

    ! Calc S_box1
    function calc_shelf_Sbox_1(to,so,T_box) result(S_box)

        implicit none

        real(prec), intent(IN)     :: to, so, T_box
        real(prec) :: S_box

        S_box = so - (so/(rho_ice_sw*lambda))*(to-T_box)

        return

    end function calc_shelf_Sbox_1

    ! Calc T and S for box n>1
    subroutine calc_shelf_TS_box_n(T_box, S_box, A_box, T_star, CC, a_pico, gamma_tstar)

        implicit none

        real(prec), intent(INOUT)  :: T_box, S_box
        real(prec), intent(IN)     :: A_box, T_star, CC
        real(prec), intent(IN)     :: a_pico, gamma_tstar

        ! Intern variables
        real(prec) :: g1, g2, s1
        real(prec) :: fac 

        g1 = A_box * gamma_tstar
        g2 = g1/(rho_ice_sw*lambda)
        s1 = S_box/(rho_ice_sw*lambda)        

        ! Temperature for Box i > 1
        fac = (CC + g1 - g2 * a_pico * S_box)
        if (abs(fac) .lt. 1e-6) then 
            write(*,*) "fac! ", fac 
            stop 
            fac = 1e-6
        end if 

        T_box = T_box +      g1 * T_star / fac
        S_box = S_box - s1 * g1 * T_star / fac

        return

    end subroutine calc_shelf_TS_box_n

    ! Compute overturning circulation
    function calc_overtuning(to,so,T_box,S_box,C_over,rho_star,alpha_pico,beta_pico) result(CC)

        implicit none

        real(prec), intent(IN)     :: T_box, S_box, to, so
        real(prec), intent(IN)     :: C_over, rho_star, alpha_pico, beta_pico
        real(prec) :: CC
        
        CC = C_over * rho_star * (beta_pico * (so - S_box) - alpha_pico * (to - T_box))

        return

    end function calc_overtuning
    
    ! equation 5 in the PICO paper.
    ! calculate pressure melting point from potential temperature
    function calc_theta_pm(so,a_pico,b_pico,c_pico,H_ice) result(theta_pm)

        implicit none

        real(prec), intent(IN)     :: so, H_ice
        real(prec), intent(IN)     :: a_pico, b_pico, c_pico
        real(prec) :: theta_pm

        theta_pm = a_pico*so + b_pico - c_pico*rho_ice*g*H_ice

        return

    end function calc_theta_pm

    ! equation 5 in the PICO paper.
    ! calculate pressure melting point from in-situ temperature
    function calc_T_pm(so,H_ice) result(T_pm)

        implicit none

        real(prec), intent(IN)     :: so, H_ice
        real(prec) :: T_pm
        real(prec) :: a_situ, b_situ, c_situ

        ! in-situ pressure melting point from Jenkins et al. 2010 paper
        a_situ = -0.0573         ! K/PSU
        b_situ = 0.0832 + 273.15 ! K
        c_situ = 7.53e-8         ! K/bar

        T_pm = a_situ*so + b_situ - c_situ*rho_ice*g*H_ice

        return

    end function calc_T_pm

    function calc_melt_rate_pico(T_box,pm_point,gamma_tstar) result(bmb)

        implicit none

        real(prec), intent(IN) :: T_box, pm_point
        real(prec), intent(IN) :: gamma_tstar
        real(prec) :: bmb        

        ! OJO: negativo?
        bmb = -1.0*(gamma_tstar / (lambda*rho_ice_sw))*(T_box - pm_point)

        return

    end function calc_melt_rate_pico

    subroutine apply_c_deep(bmb,c_deep,depth_deep,mask_ocn,z_bed,z_sl,n_smth)
        ! Apply c_deep for killing ice in the deep ocean 

        implicit none

        real(prec), intent(INOUT)  :: bmb(:,:)
        real(prec), intent(IN)     :: c_deep, depth_deep
        integer,    intent(IN)     :: mask_ocn(:,:)
        real(prec), intent(IN)     :: z_bed(:,:)
        real(prec), intent(IN)     :: z_sl(:,:)
        integer,    intent(IN)     :: n_smth       ! Smoothing neighborhood radius in grid points

        ! Local variables 
        integer :: i, j, nx, ny
        real(prec) :: n_pts
        logical,    allocatable :: is_c_deep(:,:)
        real(prec), allocatable :: bmb_tmp(:,:)

        nx = size(bmb,1)
        ny = size(bmb,2)

        allocate(is_c_deep(nx,ny))
        allocate(bmb_tmp(nx,ny))

        if (n_smth .lt. 0 .or. n_smth .gt. 2) then
            write(*,*) "apply_c_deep:: Error: n_smth should be 0, 1, or 2."
            write(*,*) "n_smth = ", n_smth
            stop
        end if

        ! Apply c_deep to appropriate regions or keep bmb, whichever is more negative

        where(mask_ocn .eq. 2 .and. z_sl-z_bed .ge. depth_deep) bmb = c_deep
        ! jablasco: parche
            !bmb = min(bmb,par%c_deep)


        ! Make sure bmb transitions smoothly to c_deep value from neighbors 

        if (n_smth .gt. 0) then

            ! Get size of neighborhood
            n_pts = real( (2*n_smth+1)**2, prec)

            bmb_tmp = bmb
            do i = 1+n_smth, nx-n_smth
            do j = 1+n_smth, ny-n_smth

                ! Apply a neighborhood average value, if deep ocean is in neighborhood of shallow ocean
                if (is_c_deep(i,j) .and. &
                    count(.not. is_c_deep(i-n_smth:i+n_smth,j-n_smth:j+n_smth)) .gt. 0) then

                    bmb(i,j) = sum(bmb_tmp(i-n_smth:i+n_smth,j-n_smth:j+n_smth))/n_pts

                end if

            end do
            end do

        end if

        return

    end subroutine apply_c_deep

end module pico_physics
