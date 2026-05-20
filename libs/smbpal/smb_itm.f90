module smb_itm
    ! This module contains all subroutines related to the calculation of
    ! surface mass balance. The external interface to this module is via
    ! the melt_budget() subroutine, which takes arrays from the main
    ! program as input and outputs arrays of the calculated variables.

    use nml
    use smbpal_precision 

    implicit none

    real(prec), parameter :: sec_day = 86400.0   ! [sec]
    real(prec), parameter :: rho_w   = 1.d3      ! Density of pure water [kg/m3]
    real(prec), parameter :: L_m     = 3.35e5    ! Latent heat of melting [J/kg]

    type itm_par_class
        real(prec) :: trans_a, trans_b, trans_c 
        real(prec) :: itm_c, itm_t, itm_b, itm_lat0  
        real(prec) :: H_snow_max
        real(prec) :: Pmaxfrac
        real(prec) :: H_snow_crit_desert
        real(prec) :: H_snow_crit_forest
        real(prec) :: melt_crit 
        real(prec) :: alb_ocean, alb_land, alb_forest, alb_ice 
        real(prec) :: alb_snow_dry, alb_snow_wet 

    end type 

    private 
    public :: itm_par_class, itm_par_load
    public :: calc_snowpack_budget_step 

contains

    subroutine itm_par_load(par,filename,init,group)

        type(itm_par_class)     :: par
        character(len=*), intent(IN) :: filename 
        logical, optional :: init 
        logical :: init_pars 
        character(len=*), intent(IN), optional   :: group

        ! Local variables 
        ! integer :: file_unit 
        character(len=32) :: nml_group
        
        ! Make sure we know the namelist group for the itm block
        if (present(group)) then
            nml_group = trim(group)
        else
            nml_group = "itm"         ! Default parameter blcok name
        end if

        ! ! Local parameter definitions (identical to object)
        ! real(prec) :: trans_a, trans_b, trans_c 
        ! real(prec) :: itm_c, itm_t, itm_b, itm_lat0 
        ! real(prec) :: H_snow_max
        ! real(prec) :: Pmaxfrac
        ! real(prec) :: H_snow_crit_desert
        ! real(prec) :: H_snow_crit_forest
        ! real(prec) :: melt_crit 
        ! real(prec) :: alb_ocean, alb_land, alb_forest, alb_ice 
        ! real(prec) :: alb_snow_dry, alb_snow_wet 

        init_pars = .FALSE.
        if (present(init)) init_pars = .TRUE.

        call nml_read(filename,nml_group,"trans_a",par%trans_a,init=init_pars)
        call nml_read(filename,nml_group,"trans_b",par%trans_b,init=init_pars)
        call nml_read(filename,nml_group,"trans_c",par%trans_c,init=init_pars)
        call nml_read(filename,nml_group,"itm_c",par%itm_c,init=init_pars)
        call nml_read(filename,nml_group,"itm_t",par%itm_t,init=init_pars)
        call nml_read(filename,nml_group,"itm_b",par%itm_b,init=init_pars)
        call nml_read(filename,nml_group,"itm_lat0",par%itm_lat0,init=init_pars)
        call nml_read(filename,nml_group,"H_snow_max",par%H_snow_max,init=init_pars)
        call nml_read(filename,nml_group,"Pmaxfrac",par%Pmaxfrac,init=init_pars)
        call nml_read(filename,nml_group,"H_snow_crit_desert",par%H_snow_crit_desert,init=init_pars)
        call nml_read(filename,nml_group,"H_snow_crit_forest",par%H_snow_crit_forest,init=init_pars)
        call nml_read(filename,nml_group,"melt_crit",par%melt_crit,init=init_pars)
        call nml_read(filename,nml_group,"alb_ocean",par%alb_ocean,init=init_pars)
        call nml_read(filename,nml_group,"alb_land",par%alb_land,init=init_pars)
        call nml_read(filename,nml_group,"alb_forest",par%alb_forest,init=init_pars)
        call nml_read(filename,nml_group,"alb_ice",par%alb_ice,init=init_pars)
        call nml_read(filename,nml_group,"alb_snow_dry",par%alb_snow_dry,init=init_pars)
        call nml_read(filename,nml_group,"alb_snow_wet",par%alb_snow_wet,init=init_pars)

        ! namelist /itm/ trans_a, trans_b, trans_c, itm_c, itm_t, itm_b, itm_lat0, &
        ! H_snow_max, Pmaxfrac, &
        ! H_snow_crit_desert, H_snow_crit_forest, melt_crit, alb_ocean, alb_land, &
        ! alb_forest, alb_ice, alb_snow_dry, alb_snow_wet 

                
        ! ! Store initial values in local parameter values 
        ! trans_a            = par%trans_a
        ! trans_b            = par%trans_b
        ! trans_c            = par%trans_c
        ! itm_c              = par%itm_c
        ! itm_t              = par%itm_t 
        ! itm_b              = par%itm_b 
        ! itm_lat0           = par%itm_lat0 
        ! H_snow_max         = par%H_snow_max
        ! Pmaxfrac           = par%Pmaxfrac
        ! H_snow_crit_desert = par%H_snow_crit_desert
        ! H_snow_crit_forest = par%H_snow_crit_forest
        ! melt_crit          = par%melt_crit
        ! alb_ocean          = par%alb_ocean
        ! alb_land           = par%alb_land
        ! alb_forest         = par%alb_forest
        ! alb_ice            = par%alb_ice
        ! alb_snow_dry       = par%alb_snow_dry
        ! alb_snow_wet       = par%alb_snow_wet

        ! ! Read parameters from input namelist file
        ! inquire(file=trim(filename),NUMBER=file_unit)
        ! if (file_unit .gt. 0) then 
        !     read(file_unit,nml=itm)
        ! else
        !     open(7,file=trim(filename))
        !     read(7,nml=itm)
        !     close(7)
        ! end if 

        ! ! Store local parameter values in output object
        ! par%trans_a            = trans_a
        ! par%trans_b            = trans_b
        ! par%trans_c            = trans_c
        ! par%itm_c              = itm_c
        ! par%itm_t              = itm_t 
        ! par%itm_b              = itm_b 
        ! par%itm_lat0           = itm_lat0 
        ! par%H_snow_max         = H_snow_max
        ! par%Pmaxfrac           = Pmaxfrac
        ! par%H_snow_crit_desert = H_snow_crit_desert
        ! par%H_snow_crit_forest = H_snow_crit_forest
        ! par%melt_crit          = melt_crit
        ! par%alb_ocean          = alb_ocean
        ! par%alb_land           = alb_land
        ! par%alb_forest         = alb_forest
        ! par%alb_ice            = alb_ice
        ! par%alb_snow_dry       = alb_snow_dry
        ! par%alb_snow_wet       = alb_snow_wet

        return

    end subroutine itm_par_load

   
    elemental subroutine calc_snowpack_budget_step(par,dt,lat,z_srf,H_ice,S,t2m,PDDs,pr,sf, &
                                           H_snow,alb_s,smbi,smb,melt,runoff,refrz,melt_net)
    ! Determine the total melt, accumulation and surface mass balance at a given point
    !  * Modified from rembo subroutine `melt_budget`
    !  * input in mm water equivalent
    ! Note: Definitions as in Ettema et al (2009) supplementary information,
    !     SMB  = sf   + rf   - runoff     [kg m2 / d ] == [mm / d]
    !   runoff = rain + melt - refrz      [kg m2 / d ] == [mm / d]

    implicit none

    type(itm_par_class), intent(IN)    :: par 
    real(prec),          intent(IN)    :: dt         ! Timestep [days], dt >= 1
    real(prec),          intent(IN)    :: lat
    real(prec),          intent(IN)    :: z_srf, H_ice, S, t2m, PDDs, pr, sf 
    real(prec),          intent(INOUT) :: H_snow  
    real(prec),          intent(OUT)   :: alb_s, smbi, smb, melt, runoff, refrz
    real(prec),          intent(OUT)   :: melt_net  

    ! Local variables
    real(prec) :: itm_c 
    real(prec) :: melt_pot
    real(prec) :: rf, atrans, rfac 
    real(prec) :: melted_snow, melted_ice, snow_to_ice 
    real(prec) :: refrz_rain, refrz_snow

    ! Determine rainfall rate from precip and snowfall [mm d-1]
    rf = pr - sf

    ! Determine preliminary surface and planetary albedo
    alb_s = calc_albedo_surface(par,z_srf,H_ice,H_snow,PDDs)

    ! Add additional snowfall [mm]
    H_snow  = H_snow + sf*dt

    ! Get amount of potential melt from ITM scheme
!     atrans   = calc_atmos_transmissivity(z_srf,H_ice,par%trans_a,par%trans_b,par%trans_c)
    atrans   = calc_atmos_transmissivity(z_srf,par%trans_a,par%trans_b)
    if (abs(par%itm_lat0) .lt. 90.0) then
        itm_c    = itm_c_lat(par%itm_c,par%itm_b,par%itm_lat0,lat)
    else 
        itm_c = par%itm_c
    end if
    melt_pot = calc_itm(S,t2m-273.15,alb_s,atrans,itm_c,par%itm_t)

    ! Determine how much snow and ice would be melted today [mm]
    if (melt_pot*dt .gt. H_snow) then 
      
      ! All snow is melted, the rest of energy converted to melt some ice
      ! The rest of energy will go into melting ice
      melted_snow = H_snow 
      melted_ice  = melt_pot*dt - H_snow
      
    else
      
      ! Snow melt will use all energy, none left for ice melt
      melted_snow = melt_pot*dt
      melted_ice  = 0.0
      
    end if    

    ! Total ablation
    melt   = melted_snow + melted_ice
    
    ! Remove any melted snow from the snow height budget
    H_snow = H_snow - melted_snow

    ! To avoid numerical issues with dt>1
    H_snow = max(H_snow,0.0)

    ! Adjust the albedo (accounting for actual amount of melt)
    alb_s = calc_albedo_surface(par,z_srf,H_ice,H_snow,PDDs,melt=melt/dt)

    ! Determine what fraction of the melted snow and rain will refreeze, 
    ! (Note: rf is zero if not on ice sheet or there is no snow cover)
!     if ( H_ice .gt. 0.0 .and. H_snow .gt. 0.0 ) then 
!     if ( H_snow .gt. 0.0 ) then 

    ! Modify refreezing factor based on height of snow 
    rfac = par%Pmaxfrac * sf/max(1e-3,pr)
    rfac = rfac + min(1.0,H_snow/1e3) * (1.0 - rfac) ! linear function increasing to rf=1 as H_snow increases to 1m.

    ! Determine the actual maximum amount of refreezing (limited to snowpack thickness)
    refrz_rain     = min(rf*dt*rfac,H_snow)                   ! First rain takes up refreezing capacity
    refrz_snow     = min(melted_snow*rfac,H_snow-refrz_rain)  ! Melted_snow uses remaining capacity if it needs it
    refrz          = refrz_snow + refrz_rain                  ! Total refreezing

    ! Any snow thickness used for refreezing turns to ice 
    ! Additionally determine how much new ice is made from compression of remaining snow
    snow_to_ice = refrz
    H_snow = H_snow - refrz
    if (H_snow .gt. par%H_snow_max) then 
      ! Assume excess contributes to new ice
      snow_to_ice = snow_to_ice + (H_snow - par%H_snow_max)
      H_snow = par%H_snow_max
    end if

    ! Determine net runoff
    runoff = (melted_snow-refrz_snow) + (rf*dt-refrz_rain) + melted_ice

    ! Get the global surface mass balance [mm]
    smb = (sf + rf)*dt - runoff

    ! Get the internal surface mass balance (what ice sheet model needs)
    ! These should apply at separate places, so
    ! smbi = -melted_ice, for negative mass balance
    ! smbi =     new_ice, for positive mass balance
    smbi = snow_to_ice + refrz - melted_ice

    ! Calculate net melt (refreezing minus total melt) for surface temp energy adjustment
    if (H_ice .gt. 0.0) then 
        melt_net = refrz - melt 
    else 
        melt_net = refrz - melted_snow    ! refrz is zero here, but keep it for consistency
    end if 

    ! Convert back to daily rates [mm/d]
    melt     = melt/dt 
    melt_net = melt_net/dt 
    runoff   = runoff/dt 
    refrz    = refrz/dt 
    smb      = smb/dt
    smbi     = smbi/dt 

    return

    end subroutine calc_snowpack_budget_step
    
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine : g e t _ a l b e d o
    ! Author     : Alex Robinson
    ! Purpose    : Determine the current surface and planetary albedos
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    elemental function calc_albedo_surface(par,z_srf,H_ice,H_snow,PDDs,melt) result(alb)
    
        implicit none
        
        type(itm_par_class), intent(IN) :: par 
        real(prec),   intent(IN) :: z_srf, H_ice, H_snow, PDDs 
        real(prec),   intent(IN), optional :: melt
        real(prec) :: alb 

        ! Local variables
        real(prec) :: H_snow_crit, depth, as_snow, alb_bg, melt_now
        integer :: n, k
        
        ! Determine the critical snow height based on number of PDDs,
        ! which correspond to vegetation growth
        ! 0-100 PDDs    = desert, 10mm
        ! 100-1000 PDDs = tundra (grass), linearly increasing between 10mm => 100mm
        ! 1000 PDDs +   = forest, 100mm
        if ( PDDs .le. 100.0 ) then 
            H_snow_crit = par%H_snow_crit_desert ! 10 mm
        else if (PDDs .le. 1000.0) then 
            H_snow_crit = par%H_snow_crit_desert +   &
                (par%H_snow_crit_forest-par%H_snow_crit_desert) * (PDDs-100.0)/(1000.0-100.0)
        else
            H_snow_crit = par%H_snow_crit_forest ! 100 mm
        end if
        
        ! Determine the scaled snowdepth
        depth = min( H_snow / H_snow_crit, 1.0 ) 
        
        ! Determine the amount of melt to affect abledo calculation.
        ! Default is high melt everywhere, so that jump in albedo is avoided
        ! at onset of melt season. After calculating actual melt, this is included 
        ! as an argument and albedo is recalculated.
        melt_now = par%melt_crit+1.0
        if ( present(melt) ) melt_now = melt
    
        ! Figure out what the surface albedo would be with no snow, 
        ! based on the type of ground underneath ( ice or land )
        if (z_srf .le. 0.0) then 
            alb_bg = par%alb_ocean

        else if (z_srf .gt. 0.0 .and. H_ice .eq. 0.0) then 
             alb_bg = par%alb_land*(1d3-min(PDDs,1d3))/(1d3-0d0) &
                                        + par%alb_forest*(min(PDDs,1d3))/(1d3-0d0)
        else 
            alb_bg = par%alb_ice 

        end if

        ! Determine current snow albedo: if melt gt eg, 1mm/day, then change albedo to melting snow!
        as_snow = par%alb_snow_dry
        if (melt_now .gt. par%melt_crit) as_snow = par%alb_snow_wet  
        
        ! Get current surface albedo to be used for ITM melt scheme
        ! It will either be the maximum albedo (that of dry snow: as_snow)
        ! or the wet snow albedo plus a fraction depending on height of snow
        ! minimum albedo now should be that of ice / wet snow minimum
        alb = alb_bg + depth*(as_snow-alb_bg)
    
        return
    
    end function calc_albedo_surface
    
    elemental function calc_albedo_planet(alb_s,a,b) result(alb_p)

        implicit none 

        real(prec), intent(IN) :: alb_s, a, b 
        real(prec) :: alb_p 

        alb_p = a + b*alb_s
    
        return 

    end function calc_albedo_planet

    elemental function calc_atmos_transmissivity(z_srf,a,b) result(at)

        implicit none 

        real(prec), intent(IN) :: z_srf, a, b 
        real(prec) :: at 

!         at = a + b*max(z_srf,0.d0)
        at = a + b*max(z_srf,0.d0)**0.5

        return 

    end function calc_atmos_transmissivity

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Subroutine : i t m (insolation temperature melt)
    ! Author     : Alex Robinson
    ! Purpose    : Determine the total melt, accumulation and 
    !              surface mass balance at a given point
    !              * output in mm water equivalent per day 
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    elemental function calc_itm(S,t2m,alb_s,atrans,c,t) result(melt)
        ! Determine the potential melt rate
        ! t2m = [degrees Celcius] !!!
        implicit none

        real(prec), intent(IN) :: S, t2m, alb_s, atrans, c, t 
        real(prec) :: melt

        ! Calculate potential melt [m/s]
        melt = (atrans*(1.d0 - alb_s)*S + c + t*t2m) / (rho_w*L_m)

        ! Convert: [m/s] => [mm/day], only positive melt
        melt = max( melt, 0.d0 ) * sec_day * 1d3  

        return

    end function calc_itm


    elemental function itm_c_lat(c,b,lat0,lat) result(c2D)
        ! Determine the potential melt rate
        ! t2m = [degrees Celcius] !!!
        implicit none

        real(prec), intent(IN) :: c, b, lat0, lat  
        real(prec) :: c2D

        ! Calculate the c parameter as a function of latitude 
        c2D = c + b*(lat-lat0)

        return

    end function itm_c_lat

end module smb_itm 

