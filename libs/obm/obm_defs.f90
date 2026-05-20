module obm_defs

    use ncio 
    use nml 

    implicit none 

    ! Choose the preciision of the library (sp,dp)
    integer,  parameter :: preci = kind(1.0)

    type obm_param_class 
        ! Box volume
        real(preci) :: vol_n   ! North box
        real(preci) :: vol_nd  ! North (deep) box
        real(preci) :: vol_t   ! Tropical box
        real(preci) :: vol_td  ! Tropical (deep) box
        real(preci) :: vol_s   ! South box
        real(preci) :: vol_sd  ! South (deep) box

        ! Box reference temperature
        real(preci) :: t_ref_n   
        real(preci) :: t_ref_nd  
        real(preci) :: t_ref_t   
        real(preci) :: t_ref_td  
        real(preci) :: t_ref_s   
        real(preci) :: t_ref_sd  

        ! Atmosphere temperatures
        real(preci) :: t_atm_n   
        real(preci) :: t_atm_t  
        real(preci) :: t_atm_s 

        ! Box relaxation time scale
        real(preci) :: tau_n   
        real(preci) :: tau_nd  
        real(preci) :: tau_t   
        real(preci) :: tau_td  
        real(preci) :: tau_s   
        real(preci) :: tau_sd  
        
        ! Box salinities
        real(preci) :: s_tot   ! [psu] Total mean salinity of North and South box 
        
        ! Other
        real(preci) :: alpha   ! [g/K??]
        real(preci) :: beta    ! [psu/K??]
        real(preci) :: mu      ! [?]

        ! Nautilus specific parameters
        real(preci) :: specific_heat_capacity
        real(preci) :: density_of_seawater
        real(preci) :: therm_exp_coeff
        real(preci) :: haline_exp_coeff
        real(preci) :: reference_salinity
        real(preci) :: overturning_threshold

        real(preci) :: thermal_ampl_north
        real(preci) :: thermal_ampl_tropics
        real(preci) :: thermal_ampl_south

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

        real(preci) :: m_init 
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
        real(preci) :: thetas_init
        real(preci) :: thetan_init
        real(preci) :: thetat_init
        real(preci) :: lambdan_init
        real(preci) :: lambdas_init
        real(preci) :: lambdat_init
        real(preci) :: lambdatd_init
    end type 

    type obm_class
        type(obm_param_class) :: par 

        real(preci) :: tn,tt,ts ! Box temperatures (surface)
        real(preci) :: tnd,ttd,tsd ! Box temperatures (deep)
        real(preci) :: sn,st,ss ! Box salinities (surface)
        real(preci) :: snd,std,ssd ! Box salinities (deep)
        real(preci) :: m    ! thc
        real(preci) :: fn,ft,fs ! external fluxes
        real(preci) :: phin,phit  ! atmospheric humidity fluxes between boxes
        real(preci) :: thetan,thetat,thetas ! Atmospheric temperatures
        real(preci) :: lambdan,lambdand,lambdat,lambdatd,lambdas,lambdasd   ! thermal coupling constants

    end type

    public

end module obm_defs