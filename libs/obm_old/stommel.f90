module stommel
    ! Jorge, a la Alex 

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
    public :: stommel_init 
    public :: stommel_update 

contains 

    subroutine stommel_init(stml, filename, group)

        implicit none 

        type(obm_class), intent(INOUT) :: stml
        character(len=*),  intent(IN) :: filename, group 

        ! Read parameters from namelist file
        call nml_read(filename, group, "vol_n", stml%par%vol_n)
        call nml_read(filename, group, "vol_s", stml%par%vol_s)
        call nml_read(filename, group, "t_ref_n", stml%par%t_ref_n)
        call nml_read(filename, group, "t_ref_t", stml%par%t_ref_t)
        call nml_read(filename, group, "t_ref_s", stml%par%t_ref_s)
        call nml_read(filename, group, "tau_n", stml%par%tau_n)
        call nml_read(filename, group, "tau_t", stml%par%tau_t)
        call nml_read(filename, group, "tau_s", stml%par%tau_s)
        call nml_read(filename, group, "s_tot", stml%par%s_tot)
        call nml_read(filename, group, "alpha", stml%par%alpha)
        call nml_read(filename, group, "beta", stml%par%beta)
        call nml_read(filename, group, "mu", stml%par%mu)

        ! Initialize prognostic states of the model
        call nml_read(filename, group, "tn_init", stml%tn)
        call nml_read(filename, group, "tt_init", stml%tt)
        call nml_read(filename, group, "ts_init", stml%ts)
        call nml_read(filename, group, "sn_init", stml%sn)
        call nml_read(filename, group, "sn_init", stml%ss)
        call nml_read(filename, group, "st_init", stml%st)
        call nml_read(filename, group, "m_init", stml%m)
        call nml_read(filename, group, "fn_init", stml%fn)
        call nml_read(filename, group, "fs_init", stml%fs) 

        return 

    end subroutine stommel_init 

    subroutine stommel_update(stml,dt)

        implicit none 

        type(obm_class), intent(INOUT) :: stml  
        real(prec), intent(IN)    :: dt 

        ! Local variables 
        real(prec) :: tn, tt, ts, sn, ss, st, m
        real(prec) :: fn, fs 
        real(prec) :: dtn, dts, dtt, dsn, rho_n, rho_t, delta_rho
        real(prec) :: vol_n, vol_s, t_ref_n, t_ref_t, t_ref_s, tau_n, tau_t, tau_s
        real(prec) :: s_tot, alpha, beta, mu 

        ! Populate local parameters 
        vol_n  = stml%par%vol_n 
        vol_s  = stml%par%vol_s 
        t_ref_n = stml%par%t_ref_n
        t_ref_t = stml%par%t_ref_t
        t_ref_s = stml%par%t_ref_s
        tau_n  = stml%par%tau_n 
        tau_t  = stml%par%tau_t 
        tau_s  = stml%par%tau_s 
        s_tot  = stml%par%s_tot 
        alpha  = stml%par%alpha 
        beta   = stml%par%beta  
        mu     = stml%par%mu  

        ! Polulate local variables 
        tn   = stml%tn 
        tt   = stml%tt 
        ts   = stml%ts 
        sn   = stml%sn 
        st   = stml%st 
        m  = stml%m 
        fn = stml%fn 
        fs = stml%fs 

        ! === Stommel model calculations =======================

        ! Update density and density gradient
        rho_n     = max(-alpha*tn + beta*sn,0.0_prec)
        rho_t     = max(-alpha*tt + beta*st,0.0_prec)
        delta_rho = rho_t - rho_n      

        ! Mass exchange between boxes 
        m = max(-mu*delta_rho,0.0_prec)

        ! Calculate derivatives of salinity and temperature

        dsn = (abs(m)*((st-sn)/vol_n)-fn)/dt
         
        dtn = (abs(m)*(tt-tn)/vol_n+(t_ref_n-tn)/tau_n)/dt

        dtt = (abs(m)*(tn-tt)/vol_n+(t_ref_t-tt)/tau_t)/dt

        dts = (abs(m)/vol_s+(t_ref_s-ts)/tau_s)/dt


        ! Update salinities and temperature 

        sn = max(sn+dsn,0.0_prec)
        st = max(2.0_prec*s_tot-sn,0.0_prec)

        ! Update temperature 

        tn = max(tn+dtn,0.0_prec)
        tt = max(tt+dtt,0.0_prec)
        ts = max(ts+dts,0.0_prec)

        ! === End Stommel model calculations =======================

        !write(*,*) tn, tt, ts, dsn, dtn, dtt, dts 

        ! Populate output variables 
        stml%tn   = tn 
        stml%tt   = tt 
        stml%ts   = ts 
        stml%sn   = sn 
        stml%st   = st 
        stml%m  = m 
        stml%fn = fn 
        stml%fs = fs 

        return 

    end subroutine stommel_update

end module stommel



