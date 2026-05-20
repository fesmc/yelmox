module nautilus
    ! ============================================================================================= !
    ! Nautilus Ocean Box Model
    ! ---------------------------------------------------------------------------------------------
    ! Description:
    !   * Fortran 90 code based on Zickfield and Rahmstorf (2003, Ocean Dynamics) model 
    !   * The original model was adapted for paleoclimate studies by Jorge Alvarez-Solas (2019)
    !   * Now I (Sergio Pérez-Montero) am generalizing the code for the coupling with Yelmo (2024)
    !
    ! ====================== !

    use ncio
    use nml
    use obm_defs

    implicit none

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: prec = sp
    integer,  parameter :: wp   = sp 

    real(wp), parameter :: sec_year  = 365.0*24.0*60.0*60.0   ! [s/a]
    real(wp), parameter :: pi        = 3.14159265359

    private
    public :: nautilus_init 
    public :: nautilus_update 
    
    contains

    subroutine nautilus_init(ntls, filename, group)

        implicit none 

        type(obm_class), intent(INOUT) :: ntls
        character(len=*),  intent(IN) :: filename, group 

        integer          :: n
        logical          :: use_restart
        character(len=256) :: ocn_restart

        ! Read parameters from namelist file
        call nml_read(filename, group, "use_restart", use_restart)
        call nml_read(filename, group, "restart", ocn_restart)

        call nml_read(filename, group, "specific_heat_capacity", ntls%par%specific_heat_capacity)
        call nml_read(filename, group, "density_of_seawater", ntls%par%density_of_seawater)
        call nml_read(filename, group, "therm_exp_coeff", ntls%par%therm_exp_coeff)
        call nml_read(filename, group, "haline_exp_coeff", ntls%par%haline_exp_coeff)
        call nml_read(filename, group, "reference_salinity", ntls%par%reference_salinity)
        call nml_read(filename, group, "overturning_threshold", ntls%par%overturning_threshold)
        call nml_read(filename, group, "thermal_coupling_constant", ntls%par%thermal_coupling_constant)
        call nml_read(filename, group, "emp_flow_k", ntls%par%emp_flow_k)
        call nml_read(filename, group, "vol_n", ntls%par%vol_n)
        call nml_read(filename, group, "vol_t", ntls%par%vol_t)
        call nml_read(filename, group, "vol_td", ntls%par%vol_td)
        call nml_read(filename, group, "vol_s", ntls%par%vol_s)
        call nml_read(filename, group, "depth_n", ntls%par%depth_n)
        call nml_read(filename, group, "depth_t", ntls%par%depth_t)
        call nml_read(filename, group, "depth_td", ntls%par%depth_td)
        call nml_read(filename, group, "depth_s", ntls%par%depth_s)

        call nml_read(filename, group, "thermal_ampl_north", ntls%par%thermal_ampl_north)
        call nml_read(filename, group, "thermal_ampl_tropics", ntls%par%thermal_ampl_tropics)
        call nml_read(filename, group, "thermal_ampl_south", ntls%par%thermal_ampl_south)

        call nml_read(filename, group, "nh_hydro_sensitivity", ntls%par%hn)
        call nml_read(filename, group, "sh_hydro_sensitivity", ntls%par%hs)
        call nml_read(filename, group, "nh_temp_constant", ntls%par%pnh)
        call nml_read(filename, group, "sh_temp_constant", ntls%par%psh)
        

        ! Initialize prognostic states of the model
        if (use_restart) then
            n = nc_size(ocn_restart,"time")
            call nc_read(ocn_restart,"tn",ntls%par%tn_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"tt",ntls%par%tt_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"ttd",ntls%par%ttd_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"ts",ntls%par%ts_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"sn",ntls%par%sn_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"st",ntls%par%st_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"std",ntls%par%std_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"ss",ntls%par%ss_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"fn",ntls%par%fn_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"ft",ntls%par%ft_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"fs",ntls%par%fs_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"phin",ntls%par%phin_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"phit",ntls%par%phit_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"thetan",ntls%par%thetan_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"thetat",ntls%par%thetat_init,start=[n],count=[1]) 
            call nc_read(ocn_restart,"thetas",ntls%par%thetas_init,start=[n],count=[1]) 
            call nml_read(filename, group, "lambdan_init", ntls%par%lambdan_init)
            call nml_read(filename, group, "lambdat_init", ntls%par%lambdat_init)
            call nml_read(filename, group, "lambdatd_init", ntls%par%lambdatd_init)  
            call nml_read(filename, group, "lambdas_init", ntls%par%lambdas_init)
        else
            call nml_read(filename, group, "tn_init", ntls%par%tn_init)
            call nml_read(filename, group, "tt_init", ntls%par%tt_init)
            call nml_read(filename, group, "ttd_init", ntls%par%ttd_init)
            call nml_read(filename, group, "ts_init", ntls%par%ts_init)
            call nml_read(filename, group, "sn_init", ntls%par%sn_init)
            call nml_read(filename, group, "st_init", ntls%par%st_init)
            call nml_read(filename, group, "std_init", ntls%par%std_init)
            call nml_read(filename, group, "ss_init", ntls%par%ss_init)
            call nml_read(filename, group, "fn_init", ntls%par%fn_init)
            call nml_read(filename, group, "ft_init", ntls%par%ft_init) 
            call nml_read(filename, group, "fs_init", ntls%par%fs_init) 
            call nml_read(filename, group, "phin_init", ntls%par%phin_init)
            call nml_read(filename, group, "phit_init", ntls%par%phit_init) 
            call nml_read(filename, group, "thetan_init", ntls%par%thetan_init)
            call nml_read(filename, group, "thetat_init", ntls%par%thetat_init) 
            call nml_read(filename, group, "thetas_init", ntls%par%thetas_init)
            call nml_read(filename, group, "lambdan_init", ntls%par%lambdan_init)
            call nml_read(filename, group, "lambdat_init", ntls%par%lambdat_init)
            call nml_read(filename, group, "lambdatd_init", ntls%par%lambdatd_init)  
            call nml_read(filename, group, "lambdas_init", ntls%par%lambdas_init)
        end if

        ntls%tn = ntls%par%tn_init
        ntls%tt = ntls%par%tt_init
        ntls%ttd = ntls%par%ttd_init
        ntls%ts = ntls%par%ts_init
        ntls%sn = ntls%par%sn_init
        ntls%st = ntls%par%st_init
        ntls%std = ntls%par%std_init
        ntls%ss = ntls%par%ss_init    
        ntls%fn = ntls%par%fn_init
        ntls%ft = ntls%par%tt_init
        ntls%fs = ntls%par%fs_init  
        ntls%phin = ntls%par%phin_init
        ntls%phit = ntls%par%phit_init    
        ntls%thetan = ntls%par%thetan_init
        ntls%thetat = ntls%par%thetat_init
        ntls%thetas = ntls%par%thetas_init    
        ntls%lambdan = ntls%par%lambdan_init
        ntls%lambdat = ntls%par%lambdat_init
        ntls%lambdatd = ntls%par%lambdatd_init
        ntls%lambdas = ntls%par%lambdas_init    

        return

    end subroutine nautilus_init

    subroutine nautilus_update(ntls,dt)

        implicit none 

        type(obm_class), intent(INOUT) :: ntls  
        real(prec), intent(IN)    :: dt 
        real(prec), parameter :: factor = 365*24*60*60 * 1e6

        ! Local variables
        real(prec) :: dTsdt, dTndt, dTtdt, dTtddt 
        real(prec) :: dSsdt, dSndt, dStdt, dStddt 
        real(prec) :: fn, ft, fs, phin, phit
        real(prec), parameter :: maximum_overturning=100.0, maximum_salinity=50.0, minimum_salinity=0.0

        ! Transform external fluxes (Sv) to m3a-1 (1 Sv = 1e6 m3s-1 = 1e6 * 365*24*60*60 m3yr-1)
        fn = ntls%fn * factor
        ft = ntls%ft * factor
        fs = ntls%fs * factor

        phin = ntls%phin * factor
        phit = ntls%phit * factor

        ! Update diagnostic variables
        ntls%m = max(calc_overturning(ntls%par%emp_flow_k, ntls%par%haline_exp_coeff, ntls%par%therm_exp_coeff, ntls%ts, ntls%ss, ntls%tn, ntls%sn), ntls%par%overturning_threshold*factor)
        call calc_all_lambda(ntls)

        ! Compute time derivatives, box loop: [...-> South --> Tropics --> North --> Tropics Deep -...]
        dTsdt = calc_temperature_derivative(ntls%par%vol_s, ntls%lambdas, ntls%m, -ntls%ts, ntls%ttd, ntls%thetas)               
        dTtdt = calc_temperature_derivative(ntls%par%vol_t, ntls%lambdat, ntls%m, -ntls%tt, ntls%ts, ntls%thetat)
        dTndt = calc_temperature_derivative(ntls%par%vol_n, ntls%lambdan, ntls%m, -ntls%tn, ntls%tt, ntls%thetan)                                          
        dTtddt = calc_temperature_derivative(ntls%par%vol_td, 0.0, ntls%m, -ntls%ttd, ntls%tn, 0.0)   
        
        dSsdt = calc_salinity_derivative(ntls%par%vol_s, ntls%par%reference_salinity, ntls%m, -ntls%ss, ntls%std, phit, -fs)
        dStdt = calc_salinity_derivative(ntls%par%vol_t, ntls%par%reference_salinity, ntls%m, -ntls%st, ntls%ss, phin, -phit)
        dSndt = calc_salinity_derivative(ntls%par%vol_n, ntls%par%reference_salinity, ntls%m, -ntls%sn, ntls%st, -fn, -phin)
        dStddt = calc_salinity_derivative(ntls%par%vol_td, 0.0, ntls%m, -ntls%std, ntls%sn, 0.0, 0.0)

        ! Update prognostic variables
        ntls%ts = update_state(ntls%ts, dTsdt, dt)
        ntls%tt = update_state(ntls%tt, dTtdt, dt)
        ntls%tn = update_state(ntls%tn, dTndt, dt)
        ntls%ttd = update_state(ntls%ttd, dTtddt, dt)

        ntls%ss = max(min(update_state(ntls%ss, dSsdt, dt), maximum_salinity), minimum_salinity)
        ntls%st = max(min(update_state(ntls%st, dStdt, dt), maximum_salinity), minimum_salinity)
        ntls%sn = max(min(update_state(ntls%sn, dSndt, dt), maximum_salinity), minimum_salinity)
        ntls%std = max(min(update_state(ntls%std, dStddt, dt), maximum_salinity), minimum_salinity)

        return

    end subroutine nautilus_update

    subroutine calc_all_lambda(ntls)
        implicit none 

        type(obm_class), intent(INOUT) :: ntls 

        ntls%lambdas = calc_individual_thermal_coupling_constant(ntls%par%thermal_coupling_constant, &
                                                                 ntls%par%specific_heat_capacity, ntls%par%density_of_seawater, &
                                                                 ntls%par%depth_s)
        ntls%lambdat = calc_individual_thermal_coupling_constant(ntls%par%thermal_coupling_constant, &
                                                                 ntls%par%specific_heat_capacity, ntls%par%density_of_seawater, &
                                                                 ntls%par%depth_t)
        ntls%lambdan = calc_individual_thermal_coupling_constant(ntls%par%thermal_coupling_constant, &
                                                                 ntls%par%specific_heat_capacity, ntls%par%density_of_seawater, &
                                                                 ntls%par%depth_n)                                  

        return

    end subroutine calc_all_lambda

    function update_state(u, dudt, deltat) result(newu)
        real(prec) :: u, dudt, deltat, newu
        newu = u + dudt * deltat
        return 
    end function update_state

    function calc_overturning(k, beta, alpha, t1, s1, t2, s2) result(m)
        real(prec) :: k, beta, alpha, t1, s1, t2, s2, m
        m = k * (beta * (s2 - s1) - alpha * (t2 - t1))
        return
    end function calc_overturning

    function calc_temperature_derivative(vol1, lambda1, m, t1, t2, theta1) result(dTdt)
        ! WARNING: make sure you provide the correct sign to the flux
        real(prec) :: vol1, lambda1, m, t1, t2, theta1, dTdt
        dTdt = m / vol1 * (t2 + t1) + lambda1 * (theta1 + t1)
        return
    end function calc_temperature_derivative

    function calc_salinity_derivative(vol1, s0, m, s1, s2, f1, f2) result(dSdt)
        ! WARNING: make sure you provide the correct sign to the flux
        real(prec) :: vol1, s0, f1, f2, m, s1, s2, dSdt
        dSdt = m / vol1 * (s2 + s1) + s0 / vol1 * (f1 + f2)
        return
    end function calc_salinity_derivative

    function calc_individual_thermal_coupling_constant(l, c, rho, z) result(lambda)
        real(preci) :: l, c, rho, z, lambda
        lambda = l / (c * rho * z)
        return
    end function calc_individual_thermal_coupling_constant

end module nautilus