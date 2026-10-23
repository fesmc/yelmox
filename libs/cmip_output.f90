module cmip_output
    ! CMIP/ISMIP-formatted output of the ice sheet ([output] write_cmip, every
    ! dt_cmip): 2D fields (yelmo_cmip.nc) and 1D integrals (yelmo_ts_cmip.nc),
    ! from Yelmo and the marine shelf. The marine-shelf fields are read on the
    ! Yelmo grid (grid_shelf == grid_ice).
    !
    !   cmip_write_init      create both files
    !   write_step_2D_cmip   CMIP/ISMIP-formatted 2D output
    !   write_step_1D_cmip   CMIP/ISMIP-formatted 1D output
    !
    ! The mass fluxes are built once (cmip_mass_fluxes) from the fields Yelmo
    ! applies; the 2D file writes them, the 1D file their integrals over the
    ! cells (times the cell area), so both files agree.

    use ncio
    use yelmo
    use marine_shelf

    implicit none

    private
    public :: cmip_write_init
    public :: write_step_2D_cmip
    public :: write_step_1D_cmip

contains

    subroutine cmip_write_init(ylmo, file2D, file1D, time)
        ! Create the 2D and 1D CMIP files; the 1D integrals cover the cells
        ! where ice is allowed.
        type(yelmo_class), intent(in) :: ylmo
        character(len=*),  intent(in) :: file2D, file1D
        real(dp),          intent(in) :: time

        call yelmo_write_init(ylmo, file2D, time_init=time, units="years")
        call yelmo_write_reg_init(ylmo, file1D, time_init=time, units="years", &
                                  mask=(ylmo%bnd%mask_ice /= MASK_ICE_NONE))
    end subroutine cmip_write_init

    subroutine write_step_2D_cmip(ylmo, mshlf, filename, time)
        ! Writes all mandatory (and key optional) 2-D ISMIP7 variables.
        ! ST = snapshot (end-of-year); FL = yearly-average flux.
        ! -------------------------------------------------------------------------
        
        implicit none
        
        type(yelmo_class),    intent(IN) :: ylmo
        type(marshelf_class), intent(IN) :: mshlf
        character(len=*),     intent(IN) :: filename
        real(dp),             intent(IN) :: time
        
        ! ---- local variables ------------------------------------------------
        integer  :: ncid, n, i, j, k, nz
        
        real(wp) :: yr_to_sec
        
        ! 2-D working arrays
        real(wp), allocatable :: acabf(:,:), libmassbfgr(:,:), libmassbffl(:,:)
        real(wp), allocatable :: licalvf(:,:), ligroundf(:,:)
        real(wp), allocatable :: z_base(:,:)
        real(wp), allocatable :: T_top_ice(:,:), T_base_grnd(:,:), T_base_flt(:,:), T_avg(:,:)
        real(wp), allocatable :: dTdz_base_grnd(:,:), dTdz_base_flt(:,:)
        real(wp), allocatable :: tfbase(:,:)
        real(wp), allocatable :: ux_aa(:,:), uy_aa(:,:) 
        real(wp), allocatable :: uz_s_masked(:,:), uz_b_masked(:,:)   
        
        ! ---- allocate -------------------------------------------------------
        allocate(z_base          (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_top_ice       (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_base_grnd     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_base_flt      (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_avg           (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(dTdz_base_grnd  (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(dTdz_base_flt   (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(tfbase          (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(ux_aa           (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(uy_aa           (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(uz_s_masked     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(uz_b_masked     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        
        ! ---- initialise -----------------------------------------------------
        z_base          = 0.0_wp;  T_top_ice        = 0.0_wp
        T_base_grnd     = 0.0_wp;  T_base_flt       = 0.0_wp
        T_avg           = 0.0_wp;  tfbase           = 0.0_wp
        dTdz_base_grnd  = 0.0_wp;  dTdz_base_flt    = 0.0_wp
        ux_aa           = 0.0_wp;  uy_aa            = 0.0_wp
        
        ! Yelmo's year [s] (velocities per second)
        yr_to_sec = ylmo%bnd%c%sec_year
        nz = ylmo%dyn%par%nz_aa
        
        ! ---- derived fields -------------------------------------------------
        
        ! Mass fluxes [kg m-2 s-1]
        call cmip_mass_fluxes(ylmo, acabf, libmassbfgr, libmassbffl, licalvf, ligroundf)
        
        ! Ice-base elevation
        z_base = ylmo%tpo%now%z_base
        
        ! Temperature fields
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp) T_top_ice = ylmo%thrm%now%T_ice(:,:,nz)
        where (ylmo%tpo%now%f_grnd .gt. 0.0_wp) T_base_grnd = ylmo%thrm%now%T_ice(:,:,1)
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp) T_base_flt = ylmo%thrm%now%T_ice(:,:,1)
        
        ! Depth-averaged temperature
        if (nz .gt. 0) then
            do k = 1, nz
                where (ylmo%tpo%now%H_ice .gt. 0.0_wp)
                    T_avg = T_avg + ylmo%thrm%now%T_ice(:,:,k)
                end where
            end do
            where (ylmo%tpo%now%H_ice .gt. 0.0_wp)
                T_avg = T_avg / real(nz, wp)
            end where
        end if
        
        ! Vertical basal temperature gradient (first-order upward difference)
        if (nz .gt. 1) then
            where (ylmo%tpo%now%f_grnd .gt. 0.0_wp .and. ylmo%tpo%now%H_ice .gt. 1.0_wp)
                dTdz_base_grnd = (ylmo%thrm%now%T_ice(:,:,2) - ylmo%thrm%now%T_ice(:,:,1)) &
                                 / (ylmo%tpo%now%H_ice / real(nz - 1, wp))
            end where
            where (ylmo%tpo%now%H_ice .gt. 1.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp)
                dTdz_base_flt  = (ylmo%thrm%now%T_ice(:,:,2) - ylmo%thrm%now%T_ice(:,:,1)) &
                                 / (ylmo%tpo%now%H_ice / real(nz - 1, wp))
            end where
        end if
        
        ! Vertical velocities (weird shape with regions)
        uz_s_masked = 0.0_wp
        uz_b_masked = 0.0_wp
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp)
            uz_s_masked = ylmo%dyn%now%uz_s / yr_to_sec
            uz_b_masked = ylmo%dyn%now%uz_b / yr_to_sec
        end where

        ! Thermal forcing at ice base (floating only)
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp)
            tfbase = mshlf%now%tf_shlf
        end where
        
        ! Mean velocities interpolated onto aa-nodes (staggered → centred)
        do j = 2, ylmo%grd%G%ny - 1
        do i = 2, ylmo%grd%G%nx - 1
            ux_aa(i,j) = 0.5_wp * (ylmo%dyn%now%ux_bar(i,j) + ylmo%dyn%now%ux_bar(i-1,j))
            uy_aa(i,j) = 0.5_wp * (ylmo%dyn%now%uy_bar(i,j) + ylmo%dyn%now%uy_bar(i,j-1))
        end do
        end do
        
        ! ---- open file & find time index ------------------------------------
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        
        ! ====================================================================
        ! 2-D ST variables  (snapshot)
        ! ====================================================================
        
        call nc_write(filename, "lithk", ylmo%tpo%now%H_ice, &
            units="m", long_name="Ice thickness", &
            standard_name="land_ice_thickness", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp)
        
        call nc_write(filename, "orog", ylmo%tpo%now%z_srf, &
            units="m", long_name="Surface elevation", &
            standard_name="surface_altitude", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp)
        
        call nc_write(filename, "topg", ylmo%bnd%z_bed, &
            units="m", long_name="Bedrock elevation", &
            standard_name="bedrock_altitude", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp)
        
        call nc_write(filename, "base", z_base, &
            units="m", long_name="Ice base elevation", &
            standard_name="base_altitude", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp)
        
        call nc_write(filename, "xvelsurf", ylmo%dyn%now%ux_s / yr_to_sec, &
            units="m s-1", long_name="Surface velocity in x", &
            standard_name="land_ice_surface_x_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp/yr_to_sec, prec_rel=1e-4_wp)
        
        call nc_write(filename, "yvelsurf", ylmo%dyn%now%uy_s / yr_to_sec, &
            units="m s-1", long_name="Surface velocity in y", &
            standard_name="land_ice_surface_y_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp/yr_to_sec, prec_rel=1e-4_wp)
        
        call nc_write(filename, "zvelsurf", uz_s_masked, &
            units="m s-1", long_name="Surface velocity in z", &
            standard_name="land_ice_surface_upward_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp/yr_to_sec, prec_rel=1e-4_wp)
        
        call nc_write(filename, "xvelbase", ylmo%dyn%now%ux_b / yr_to_sec, &
            units="m s-1", long_name="Basal velocity in x", &
            standard_name="land_ice_basal_x_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp/yr_to_sec, prec_rel=1e-4_wp)
        
        call nc_write(filename, "yvelbase", ylmo%dyn%now%uy_b / yr_to_sec, &
            units="m s-1", long_name="Basal velocity in y", &
            standard_name="land_ice_basal_y_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp/yr_to_sec, prec_rel=1e-4_wp)
        
        call nc_write(filename, "zvelbase", uz_b_masked, &
            units="m s-1", long_name="Basal velocity in z", &
            standard_name="land_ice_basal_upward_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp/yr_to_sec, prec_rel=1e-4_wp)
        
        call nc_write(filename, "xvelmean", ux_aa / yr_to_sec, &
            units="m s-1", long_name="Mean velocity in x", &
            standard_name="land_ice_vertical_mean_x_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp/yr_to_sec, prec_rel=1e-4_wp)
        
        call nc_write(filename, "yvelmean", uy_aa / yr_to_sec, &
            units="m s-1", long_name="Mean velocity in y", &
            standard_name="land_ice_vertical_mean_y_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp/yr_to_sec, prec_rel=1e-4_wp)
        
        call nc_write(filename, "litemptop", T_top_ice, &
            units="K", long_name="Surface temperature", &
            standard_name="temperature_at_top_of_ice_sheet_model", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp)
        
        call nc_write(filename, "litempavg", T_avg, &
            units="K", long_name="Depth-averaged ice temperature", &
            standard_name="land_ice_temperature", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp)
        
        call nc_write(filename, "litempbotgr", T_base_grnd, &
            units="K", long_name="Basal temperature beneath grounded ice sheet", &
            standard_name="temperature_at_base_of_ice_sheet_model", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp)
        
        call nc_write(filename, "litempbotfl", T_base_flt, &
            units="K", long_name="Basal temperature beneath floating ice shelf", &
            standard_name="temperature_at_base_of_ice_sheet_model", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp)
        
        call nc_write(filename, "litempgradgr", dTdz_base_grnd, &
            units="K m-1", long_name="Vertical basal temperature gradient beneath grounded ice sheet", &
            standard_name="", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_rel=1e-4_wp)
        
        call nc_write(filename, "litempgradfl", dTdz_base_flt, &
            units="K m-1", long_name="Vertical basal temperature gradient beneath floating ice shelf", &
            standard_name="", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_rel=1e-4_wp)
        
        call nc_write(filename, "strbasemag", ylmo%dyn%now%taub, &
            units="Pa", long_name="Basal drag", &
            standard_name="land_ice_basal_drag", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_rel=1e-4_wp)
        
        call nc_write(filename, "sftgif", ylmo%tpo%now%f_ice, &
            units="1", long_name="Land ice area fraction", &
            standard_name="land_ice_area_fraction", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-4_wp)
        
        call nc_write(filename, "sftgrf", ylmo%tpo%now%f_grnd, &
            units="1", long_name="Grounded ice sheet area fraction", &
            standard_name="grounded_ice_sheet_area_fraction", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-4_wp)
        
        call nc_write(filename, "sftflf", MAX(ylmo%tpo%now%f_ice - ylmo%tpo%now%f_grnd, 0.0_wp), &
            units="1", long_name="Floating ice sheet area fraction", &
            standard_name="floating_ice_shelf_area_fraction", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-4_wp)
        
        ! ====================================================================
        ! 2-D FL variables  (yearly-average flux)
        ! ====================================================================
        
        ! Time-invariant in most setups; written without time dimension.
        ! Yelmo's Q_geo is in [mW m-2].
        call nc_write(filename, "hfgeoubed", ylmo%bnd%Q_geo * 1.0e-3_wp, &
            units="W m-2", long_name="Geothermal heat flux", &
            standard_name="upward_geothermal_heat_flux_in_land_ice", &
            dim1="xc", dim2="yc", start=[1,1], ncid=ncid)
        
        call nc_write(filename, "acabf", acabf, &
            units="kg m-2 s-1", long_name="Surface mass balance flux", &
            standard_name="land_ice_surface_specific_mass_balance_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_rel=1e-4_wp)
        
        call nc_write(filename, "libmassbfgr", libmassbfgr, &
            units="kg m-2 s-1", long_name="Basal mass balance flux beneath grounded ice", &
            standard_name="land_ice_basal_specific_mass_balance_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_rel=1e-4_wp)
        
        call nc_write(filename, "libmassbffl", libmassbffl, &
            units="kg m-2 s-1", long_name="Basal mass balance flux beneath floating ice", &
            standard_name="land_ice_basal_specific_mass_balance_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_rel=1e-4_wp)
        
        call nc_write(filename, "dlithkdt", ylmo%tpo%now%dHidt / yr_to_sec, &
            units="m s-1", long_name="Ice thickness imbalance", &
            standard_name="tendency_of_land_ice_thickness", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-4_wp/yr_to_sec)
        
        ! Calving mass balance applied by Yelmo (cmb_flt + cmb_grnd)
        call nc_write(filename, "licalvf", licalvf, &
            units="kg m-2 s-1", long_name="Calving flux", &
            standard_name="land_ice_specific_mass_flux_due_to_calving", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_rel=1e-4_wp)
        
        ! ligroundf : Grounding-line flux                          [MANDATORY]
        call nc_write(filename, "ligroundf", ligroundf, &
            units="kg m-2 s-1", long_name="Grounding line flux", &
            standard_name="land_ice_specific_grounding_line_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_rel=1e-4_wp)
        
        ! tfbase : Thermal forcing at ice base, floating           [optional]
        call nc_write(filename, "tfbase", tfbase, &
            units="K", long_name="Thermal forcing at the ice base", &
            standard_name="", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid, prec_abs=1e-3_wp)
        
        ! ---- close ----------------------------------------------------------
        call nc_close(ncid)
        
        return
        
    end subroutine write_step_2D_cmip
        
        
    subroutine write_step_1D_cmip(ylmo, mshlf, filename, time)
        ! Writes all mandatory scalar ISMIP7 (ISM_2026) variables. The fluxes
        ! are the integrals of the 2D fluxes (cmip_mass_fluxes) times the cell
        ! area: the surface mass balance over the ice (H_ice > 0), the others
        ! over the cells where they are defined.
        ! -------------------------------------------------------------------------
        
        implicit none
        
        type(yelmo_class),    intent(IN) :: ylmo
        type(marshelf_class), intent(IN) :: mshlf
        character(len=*),     intent(IN) :: filename
        real(dp),             intent(IN) :: time
        
        ! ---- local variables ------------------------------------------------
        type(yregions_class) :: reg
        
        integer  :: ncid, n
        real(wp) :: rho_ice
        real(wp) :: area             ! cell area [m2]
        
        ! 2-D mass fluxes [kg m-2 s-1]
        real(wp), allocatable :: acabf(:,:), libmassbfgr(:,:), libmassbffl(:,:)
        real(wp), allocatable :: licalvf(:,:), ligroundf(:,:)
        
        rho_ice = ylmo%bnd%c%rho_ice
        area    = ylmo%grd%G%dx * ylmo%grd%G%dy
        
        ! ---- regional object (pre-computed by Yelmo) ------------------------
        reg = ylmo%reg
        
        ! ---- the mass fluxes, integrated below -------------------------------
        call cmip_mass_fluxes(ylmo, acabf, libmassbfgr, libmassbffl, licalvf, ligroundf)
        
        ! ---- open file & find time index ------------------------------------
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        
        ! ====================================================================
        ! Scalar ST variables  (snapshot)
        ! ====================================================================
        
        ! lim : Total ice mass                                     [MANDATORY]
        call nc_write(filename, "lim", reg%V_ice * rho_ice * 1.0e9_wp, &
            units="kg", long_name="Total ice mass", &
            standard_name="land_ice_mass", &
            dim1="time", start=[n], ncid=ncid)
        
        ! limnsw : Mass above floatation                           [MANDATORY]
        call nc_write(filename, "limnsw", reg%V_sl * rho_ice * 1.0e9_wp, &
            units="kg", long_name="Mass above floatation", &
            standard_name="land_ice_mass_not_displacing_sea_water", &
            dim1="time", start=[n], ncid=ncid)
        
        ! iareagr : Grounded ice area                              [MANDATORY]
        call nc_write(filename, "iareagr", reg%A_ice_g * 1.0e6_wp, &
            units="m2", long_name="Grounded ice area", &
            standard_name="grounded_ice_sheet_area", &
            dim1="time", start=[n], ncid=ncid)
        
        ! iareafl : Floating ice area                              [MANDATORY]
        call nc_write(filename, "iareafl", reg%A_ice_f * 1.0e6_wp, &
            units="m2", long_name="Floating ice area", &
            standard_name="floating_ice_shelf_area", &
            dim1="time", start=[n], ncid=ncid)
        
        ! ====================================================================
        ! Scalar FL variables  (yearly-average flux)
        ! ====================================================================
        
        ! tendacabf : Total SMB flux                               [MANDATORY]
        call nc_write(filename, "tendacabf", sum(acabf, mask=ylmo%tpo%now%H_ice .gt. 0.0_wp) * area, &
            units="kg s-1", long_name="Total SMB flux", &
            standard_name="tendency_of_land_ice_mass_due_to_surface_mass_balance", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendlibmassbfgr : Total BMB flux, grounded               [MANDATORY]
        call nc_write(filename, "tendlibmassbfgr", sum(libmassbfgr) * area, &
            units="kg s-1", long_name="Total BMB flux beneath grounded ice", &
            standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendlibmassbffl : Total BMB flux, floating               [MANDATORY]
        call nc_write(filename, "tendlibmassbffl", sum(libmassbffl) * area, &
            units="kg s-1", long_name="Total BMB flux beneath floating ice", &
            standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendlicalvf : Total calving flux                         [MANDATORY]
        call nc_write(filename, "tendlicalvf", sum(licalvf) * area, &
            units="kg s-1", long_name="Total calving flux", &
            standard_name="tendency_of_land_ice_mass_due_to_calving", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendlifmassbf : Total ice-front melt flux                [MANDATORY]
        ! ISMIP7-2026 separates this from calving. The kinematic approach does
        ! TODO: replace with sum(ylmo%tpo%now%fmb * ...) when field is confirmed.
        call nc_write(filename, "tendlifmassbf", 0.0_wp, &
            units="kg s-1", long_name="Total ice front melting flux", &
            standard_name="tendency_of_land_ice_mass_due_to_ice_front_melting", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendligroundf : Total grounding-line flux                [MANDATORY]
        call nc_write(filename, "tendligroundf", sum(ligroundf) * area, &
            units="kg s-1", long_name="Total grounding line flux", &
            standard_name="tendency_of_grounded_ice_sheet_mass", &
            dim1="time", start=[n], ncid=ncid)
        
        ! ---- close ----------------------------------------------------------
        call nc_close(ncid)
        
        return
        
    end subroutine write_step_1D_cmip

    subroutine cmip_mass_fluxes(ylmo, acabf, libmassbfgr, libmassbffl, licalvf, ligroundf)
        ! The 2D mass fluxes of the CMIP output [kg m-2 s-1], from the fields
        ! Yelmo applies (rates of ice [m yr-1] times rho_ice/sec_year, Yelmo's
        ! own constants), the same for the 2D file and the 1D integrals:
        !   acabf        surface mass balance (tpo smb), everywhere
        !   libmassbfgr  basal mass balance (tpo bmb) of the grounded ice
        !                (H_ice > 0, f_grnd > 0), 0 elsewhere
        !   libmassbffl  basal mass balance (tpo bmb) of the floating ice
        !                (H_ice > 0, f_grnd = 0), 0 elsewhere
        !   licalvf      calving mass balance (cmb_flt + cmb_grnd)
        !   ligroundf    grounding-line flux: the kinematic flux uxy_bar*H_ice*rho_ice
        !                [kg m-1 s-1] through a cell of width dx, per cell area, in
        !                the grounded cells with mask_grz = 0 (as yelmo_calving),
        !                0 elsewhere
        type(yelmo_class),     intent(IN)  :: ylmo
        real(wp), allocatable, intent(OUT) :: acabf(:,:), libmassbfgr(:,:), libmassbffl(:,:)
        real(wp), allocatable, intent(OUT) :: licalvf(:,:), ligroundf(:,:)

        real(wp) :: conv             ! [kg m-2 s-1] per [m yr-1] of ice
        real(wp) :: dx
        logical, allocatable :: mask_grnd(:,:), mask_flt(:,:), mask_grl(:,:)

        conv = ylmo%bnd%c%rho_ice / ylmo%bnd%c%sec_year
        dx   = ylmo%grd%G%dx

        associate(tpo => ylmo%tpo%now)

        mask_grnd = (tpo%H_ice .gt. 0.0_wp .and. tpo%f_grnd .gt. 0.0_wp)
        mask_flt  = (tpo%H_ice .gt. 0.0_wp .and. tpo%f_grnd .eq. 0.0_wp)
        mask_grl  = (mask_grnd .and. tpo%mask_grz .eq. 0)

        acabf       = tpo%smb * conv
        libmassbfgr = merge(tpo%bmb, 0.0_wp, mask_grnd) * conv
        libmassbffl = merge(tpo%bmb, 0.0_wp, mask_flt)  * conv
        licalvf     = (tpo%cmb_flt + tpo%cmb_grnd) * conv
        ligroundf   = merge(ylmo%dyn%now%uxy_bar * tpo%H_ice, 0.0_wp, mask_grl) * conv / dx

        end associate

    end subroutine cmip_mass_fluxes

end module cmip_output
