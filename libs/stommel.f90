module stommel
    ! Jorge, a la Alex 

    use yelmo_defs, only : sp, dp, prec, pi 

    implicit none

    type stommel_param_class 
        real(prec) :: vol_n   ! Volume of the ocean in the North box
        real(prec) :: vol_s   ! Volume of the ocean in the South box
        real(prec) :: tn_ref  ! Reference temp. in the North  
        real(prec) :: tt_ref  ! Reference temp. in the Tropics  
        real(prec) :: ts_ref  ! Reference temp. in the South  
        real(prec) :: tau_n   ! [a] Relaxation time scale of the North box  
        real(prec) :: tau_t   ! [a] Relaxation time scale of the Tropics box 
        real(prec) :: tau_s   ! [a] Relaxation time scale of the South box 
        real(prec) :: s_tot   ! [psu] Total mean salinity of North and South box 
        real(prec) :: alpha   ! [g/K??]
        real(prec) :: beta    ! [psu/K??]
        real(prec) :: mu      ! [?]
    end type 

    type stommel_class 

        type(stommel_param_class) :: par 

        real(prec) :: tn,ts,tt
        real(prec) :: sn,ss,st
        real(prec) :: thc,fwfn,fwfs

    end type

    private
    public :: stommel_class
    public :: stommel_init 
    public :: calc_stommel 

contains 

    subroutine stommel_init(stm)

        implicit none 

        type(stommel_class), intent(INOUT) :: stm 

        ! Initialize parameters 
        stm%par%vol_n   =  300.0    ! Volume of the ocean in the North box
        stm%par%vol_s   = 3000.0    ! Volume of the ocean in the South box
        stm%par%tn_ref  =    5.0    ! Reference temp. in the North  
        stm%par%tt_ref  =   15.0    ! Reference temp. in the Tropics  
        stm%par%ts_ref  =    5.0    ! Reference temp. in the South  
        stm%par%tau_n   =   50.0    ! [a] Relaxation time scale of the North box  
        stm%par%tau_t   =  100.0    ! [a] Relaxation time scale of the Tropics box 
        stm%par%tau_s   = 1000.0    ! [a] Relaxation time scale of the South box 
        stm%par%s_tot   =   45.0    ! [psu] Total mean salinity of North and South box 
        stm%par%alpha   =  2.5e-4   ! [g/K??]
        stm%par%beta    =  0.8e-3   ! [psu/K??]
        stm%par%mu      =  1.8e4    ! [?]

        ! Initialize variables 
        stm%tn          = 0.0 
        stm%tt          = 0.0 
        stm%ts          = 0.0 
        stm%sn          = 45.0 
        stm%st          = 45.0 
        stm%thc         = 0.0 
        stm%fwfn        = 0.0 
        stm%fwfs        = 0.0 

        return 

    end subroutine stommel_init 

    subroutine calc_stommel(stm,dt)

        implicit none 

        type(stommel_class), intent(INOUT) :: stm  
        real(prec), intent(IN)    :: dt 

        ! Local variables 
        real(prec) :: tn, tt, ts, sn, st, thc
        real(prec) :: fwfn, fwfs 
        real(prec) :: dtn, dts, dtt, dsn, rho_n, rho_t, delta_rho
        real(prec) :: vol_n, vol_s, tn_ref, tt_ref, ts_ref, tau_n, tau_t, tau_s
        real(prec) :: s_tot, alpha, beta, mu 

        ! Populate local parameters 
        vol_n  = stm%par%vol_n 
        vol_s  = stm%par%vol_s 
        tn_ref = stm%par%tn_ref
        tt_ref = stm%par%tt_ref
        ts_ref = stm%par%ts_ref
        tau_n  = stm%par%tau_n 
        tau_t  = stm%par%tau_t 
        tau_s  = stm%par%tau_s 
        s_tot  = stm%par%s_tot 
        alpha  = stm%par%alpha 
        beta   = stm%par%beta  
        mu     = stm%par%mu  

        ! Polulate local variables 
        tn   = stm%tn 
        tt   = stm%tt 
        ts   = stm%ts 
        sn   = stm%sn 
        st   = stm%st 
        thc  = stm%thc 
        fwfn = stm%fwfn 
        fwfs = stm%fwfs 

        ! === Stommel model calculations =======================

        ! Update density and density gradient
        rho_n     = max(-alpha*tn + beta*sn,0.0_prec)
        rho_t     = max(-alpha*tt + beta*st,0.0_prec)
        delta_rho = rho_t - rho_n      

        ! Mass exchange between boxes 
        thc = max(-mu*delta_rho,0.0_prec)

        ! Calculate derivatives of salinity and temperature

        dsn = (abs(thc)*((st-sn)/vol_n)-fwfn)/dt
         
        dtn = (abs(thc)*(tt-tn)/vol_n+(tn_ref-tn)/tau_n)/dt

        dtt = (abs(thc)*(tn-tt)/vol_n+(tt_ref-tt)/tau_t)/dt

        dts = (abs(thc)/vol_s+(ts_ref-ts)/tau_s)/dt


        ! Update salinities and temperature 

        sn = max(sn+dsn,0.0_prec)
        st = max(2.0_prec*s_tot-sn,0.0_prec)

        ! Update temperature 

        tn = max(tn+dtn,0.0_prec)
        tt = max(tt+dtt,0.0_prec)
        ts = max(ts+dts,0.0_prec)

        ! === End Stommel model calculations =======================

        !write(*,*) tn, tt, ts, dsn, dtn, dtt, dts 

        ! Polulate output variables 
        stm%tn   = tn 
        stm%tt   = tt 
        stm%ts   = ts 
        stm%sn   = sn 
        stm%st   = st 
        stm%thc  = thc 
        stm%fwfn = fwfn 
        stm%fwfs = fwfs 

        return 

    end subroutine calc_stommel

end module stommel



