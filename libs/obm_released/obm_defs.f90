module obm_defs

    use ncio 
    use nml 

    implicit none 

    ! Choose the preciision of the library (sp,dp)
    integer,  parameter :: preci = kind(1.0)

    type obm_param_class 
        ! Box volume
        real(preci) :: vol_n   ! North box
        real(preci) :: vol_t   ! Tropical box
        real(preci) :: vol_td  ! Tropical (deep) box
        real(preci) :: vol_s   ! South box

        ! Nautilus parameters
        real(preci) :: specific_heat_capacity
        real(preci) :: density_of_seawater
        real(preci) :: therm_exp_coeff
        real(preci) :: haline_exp_coeff
        real(preci) :: reference_salinity
        real(preci) :: overturning_threshold

        real(preci) :: hn
        real(preci) :: hs
        real(preci) :: pnh
        real(preci) :: psh

        real(preci) :: thermal_coupling_constant
        real(preci) :: emp_flow_k
        real(preci) :: depth_s
        real(preci) :: depth_n 
        real(preci) :: depth_t
        real(preci) :: depth_td

        real(preci) :: tn_init
        real(preci) :: tt_init
        real(preci) :: ttd_init
        real(preci) :: ts_init
        real(preci) :: sn_init
        real(preci) :: st_init
        real(preci) :: std_init
        real(preci) :: ss_init
        real(preci) :: fn_init
        real(preci) :: ft_init
        real(preci) :: fs_init
        real(preci) :: phit_init
        real(preci) :: phin_init
        real(preci) :: tstars_init
        real(preci) :: tstarn_init
        real(preci) :: tstart_init
        real(preci) :: lambdan_init
        real(preci) :: lambdas_init
        real(preci) :: lambdat_init
        real(preci) :: lambdatd_init
    end type 

    type obm_class
        type(obm_param_class) :: par 

        real(preci) :: tn,tt,ts ! Box temperatures (surface) [degC]
        real(preci) :: ttd      ! Box temperature (tropical deep) [degC]
        real(preci) :: sn,st,ss ! Box salinities (surface) [psu]
        real(preci) :: std      ! Box salinity (tropical deep) [psu]
        real(preci) :: m        ! thc [m3/yr]
        real(preci) :: fn,ft,fs ! external fluxes [Sv]
        real(preci) :: phin,phit  ! atmospheric humidity fluxes between boxes [Sv]
        real(preci) :: tstarn,tstart,tstars ! Atmospheric temperatures [degC]
        real(preci) :: lambdan,lambdat,lambdatd,lambdas   ! thermal coupling constants

    end type

    public

end module obm_defs