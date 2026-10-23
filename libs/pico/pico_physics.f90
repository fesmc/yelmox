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

    ! Physical constants are no longer fixed here: pico_init puts them in
    ! pico_param_class from the phys_const_class it is given, and the routines
    ! below take what they need as arguments. Two groupings appear:
    !   rho_ice_g  = rho_ice*g            [Pa/m] ice overburden gradient
    !   lambda_rho = (L_ice/cp_ocn)*rho_ice/rho_sw

    private
    public :: calc_Tstar
    public :: calc_shelf_Tbox_1
    public :: calc_shelf_Sbox_1
    public :: calc_shelf_TS_box_n
    public :: calc_overtuning
    public :: calc_theta_pm
    public :: calc_melt_rate_pico

contains 

    ! =======================================
    ! 
    ! Subroutines to compute pico physics 
    ! 
    ! =======================================

    function calc_Tstar(to,so,H_ice,a_pico,b_pico,c_pico,rho_ice_g) result(T_star)

        implicit none

        real(prec), intent(IN)  :: to, so, H_ice
        real(prec), intent(IN)  :: a_pico, b_pico, c_pico
        real(prec), intent(IN)  :: rho_ice_g    ! [Pa/m] rho_ice*g
        real(prec) :: T_star

        T_star = a_pico*so + b_pico - c_pico*rho_ice_g*H_ice - to
        ! Ensures that temperature input for grounding line box should not be below pressure melting point
        ! This ensures that later equations are well solvable.
        if(T_star .gt. 0.0) T_star = 0.0 !-0.0001

        return

    end function calc_Tstar

    ! First T_box1
    function calc_shelf_Tbox_1(to,so,A_box,T_star,C_over,rho_star,gamma_tstar,alpha_pico,beta_pico, &
                              lambda_rho) result(T_box)

        implicit none

        real(prec), intent(IN)  :: A_box, T_star, to, so
        real(prec), intent(IN)  :: C_over, rho_star, gamma_tstar, alpha_pico, beta_pico
        real(prec), intent(IN)  :: lambda_rho   ! [K] (L_ice/cp_ocn)*rho_ice/rho_sw
        real(prec) :: T_box

        ! Internal variables
        real(prec) :: g1,s1,p,q,D

        g1 = gamma_tstar*A_box
        s1 = so/lambda_rho
        p = g1 / (C_over * rho_star * (beta_pico * s1 - alpha_pico))
        q = p * T_star

        D = 0.25*p*p - q
        if (D .lt. 0.0) D = 0.0
      
        T_box = to - (-0.5*p + SQRT(D))

        return

    end function calc_shelf_Tbox_1

    ! Calc S_box1
    function calc_shelf_Sbox_1(to,so,T_box,lambda_rho) result(S_box)

        implicit none

        real(prec), intent(IN)     :: to, so, T_box
        real(prec), intent(IN)     :: lambda_rho
        real(prec) :: S_box

        S_box = so - (so/lambda_rho)*(to-T_box)

        return

    end function calc_shelf_Sbox_1

    ! Calc T and S for box n>1: on input T_box, S_box are the mean
    ! temperature and salinity of box n-1, on output those of box n.
    ! Melting cools and freshens: S_n = S_n-1 - (S_n-1/lambda_rho)*(T_n-1 - T_n)
    subroutine calc_shelf_TS_box_n(T_box, S_box, A_box, T_star, CC, a_pico, gamma_tstar, lambda_rho)

        implicit none

        real(prec), intent(INOUT)  :: T_box, S_box
        real(prec), intent(IN)     :: A_box, T_star, CC
        real(prec), intent(IN)     :: a_pico, gamma_tstar
        real(prec), intent(IN)     :: lambda_rho

        ! Intern variables
        real(prec) :: g1, g2, s1
        real(prec) :: fac 

        g1 = A_box * gamma_tstar
        g2 = g1/lambda_rho
        s1 = S_box/lambda_rho        

        ! Temperature for Box i > 1
        fac = (CC + g1 - g2 * a_pico * S_box)
        if (abs(fac) .lt. 1e-6) then 
            write(*,*) "calc_shelf_TS_box_n:: Error: denominator too small, fac = ", fac 
            error stop 1
        end if 

        T_box = T_box +      g1 * T_star / fac
        S_box = S_box + s1 * g1 * T_star / fac

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
    function calc_theta_pm(so,a_pico,b_pico,c_pico,H_ice,rho_ice_g) result(theta_pm)

        implicit none

        real(prec), intent(IN)     :: so, H_ice
        real(prec), intent(IN)     :: a_pico, b_pico, c_pico
        real(prec), intent(IN)     :: rho_ice_g
        real(prec) :: theta_pm

        theta_pm = a_pico*so + b_pico - c_pico*rho_ice_g*H_ice

        return

    end function calc_theta_pm

    function calc_melt_rate_pico(T_box,pm_point,gamma_tstar,lambda_rho) result(bmb)

        implicit none

        real(prec), intent(IN) :: T_box, pm_point
        real(prec), intent(IN) :: gamma_tstar
        real(prec), intent(IN) :: lambda_rho
        real(prec) :: bmb        

        ! OJO: negativo?
        bmb = -1.0*(gamma_tstar / lambda_rho)*(T_box - pm_point)

        return

    end function calc_melt_rate_pico

end module pico_physics
