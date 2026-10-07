module cmip_output
    ! CMIP/ISMIP-formatted output of the ice sheet ([output] write_cmip, every
    ! dt_cmip): 2D fields (yelmo_cmip.nc) and 1D integrals (yelmo_ts_cmip.nc),
    ! from Yelmo and the marine shelf. The marine-shelf fields are read on the
    ! Yelmo grid (grid_shelf == grid_ice).
    !
    !   cmip_write_init      create both files
    !   write_step_2D_cmip   CMIP/ISMIP-formatted 2D output
    !   write_step_1D_cmip   CMIP/ISMIP-formatted 1D output

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
        real(wp),          intent(in) :: time

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
        real(wp),             intent(IN) :: time
        
        ! ---- local variables ------------------------------------------------
        integer  :: ncid, n, i, j, k, nz
        
        real(wp) :: rho_ice
        real(wp) :: esm_correction   ! [kg m-2 s-1] per [m yr-1]
        real(wp) :: yr_to_sec
        
        ! 2-D working arrays
        real(wp), allocatable :: bmb_grnd_masked(:,:), bmb_shlf_masked(:,:)
        real(wp), allocatable :: z_base(:,:)
        real(wp), allocatable :: T_top_ice(:,:), T_base_grnd(:,:), T_base_flt(:,:), T_avg(:,:)
        real(wp), allocatable :: dTdz_base_grnd(:,:), dTdz_base_flt(:,:)
        real(wp), allocatable :: flux_grl_2d(:,:), flux_clv_2d(:,:), tfbase(:,:)
        real(wp), allocatable :: ux_aa(:,:), uy_aa(:,:) 
        real(wp), allocatable :: uz_s_masked(:,:), uz_b_masked(:,:)   
        
        ! ---- allocate -------------------------------------------------------
        allocate(bmb_grnd_masked (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(bmb_shlf_masked (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(z_base          (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_top_ice       (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_base_grnd     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_base_flt      (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(T_avg           (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(dTdz_base_grnd  (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(dTdz_base_flt   (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(flux_grl_2d     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(flux_clv_2d     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(tfbase          (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(ux_aa           (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(uy_aa           (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(uz_s_masked     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(uz_b_masked     (ylmo%grd%G%nx, ylmo%grd%G%ny))
        
        ! ---- initialise -----------------------------------------------------
        bmb_grnd_masked = 0.0_wp;  bmb_shlf_masked = 0.0_wp
        z_base          = 0.0_wp;  T_top_ice        = 0.0_wp
        T_base_grnd     = 0.0_wp;  T_base_flt       = 0.0_wp
        T_avg           = 0.0_wp;  tfbase           = 0.0_wp
        dTdz_base_grnd  = 0.0_wp;  dTdz_base_flt    = 0.0_wp
        flux_grl_2d     = 0.0_wp;  flux_clv_2d      = 0.0_wp  
        ux_aa           = 0.0_wp;  uy_aa            = 0.0_wp
        
        ! === Unit conversion factors =========================================
        ! Taken from Yelmo's own constants, so the output cannot disagree with
        ! the model that produced it. esm_correction converts a rate in
        ! [<unit> yr-1] of ice to [kg ... s-1]; it collapses to rho_ice/sec_year,
        ! since the old two-step form was (rho_w/sec_year)*(rho_ice/rho_w) with
        ! rho_w/sec_year rounded to 3.2e-5 (1 % high) and sec_year written as
        ! 31556952 (the Gregorian year) rather than Yelmo's own.
        rho_ice        = ylmo%bnd%c%rho_ice
        yr_to_sec      = ylmo%bnd%c%sec_year
        esm_correction = rho_ice / yr_to_sec
        nz = ylmo%dyn%par%nz_aa
        
        ! ---- derived fields -------------------------------------------------
        
        ! Ice-base elevation
        z_base = ylmo%tpo%now%z_base
        
        ! BMB masked to grounded / floating
        where (ylmo%tpo%now%f_grnd .gt. 0.0_wp) bmb_grnd_masked = ylmo%thrm%now%bmb_grnd
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp) bmb_shlf_masked = ylmo%bnd%bmb_shlf
        
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

        ! Grounding-line flux (2-D, for ligroundf field)
        ! Use mask_grz convention from yelmo_calving.f90: grounded + mask_grz==0
        where (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .gt. 0.0_wp &
            .and. ylmo%tpo%now%mask_grz .eq. 0.0_wp)
            flux_grl_2d = ylmo%dyn%now%uxy_bar * ylmo%tpo%now%H_ice * rho_ice / yr_to_sec
        end where
        
        flux_clv_2d = (ylmo%tpo%now%cmb_flt+ylmo%tpo%now%cmb_grnd) * esm_correction

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
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "orog", ylmo%tpo%now%z_srf, &
            units="m", long_name="Surface elevation", &
            standard_name="surface_altitude", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "topg", ylmo%bnd%z_bed, &
            units="m", long_name="Bedrock elevation", &
            standard_name="bedrock_altitude", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "base", z_base, &
            units="m", long_name="Ice base elevation", &
            standard_name="base_altitude", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "xvelsurf", ylmo%dyn%now%ux_s / yr_to_sec, &
            units="m s-1", long_name="Surface velocity in x", &
            standard_name="land_ice_surface_x_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "yvelsurf", ylmo%dyn%now%uy_s / yr_to_sec, &
            units="m s-1", long_name="Surface velocity in y", &
            standard_name="land_ice_surface_y_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "zvelsurf", uz_s_masked, &
            units="m s-1", long_name="Surface velocity in z", &
            standard_name="land_ice_surface_upward_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "xvelbase", ylmo%dyn%now%ux_b / yr_to_sec, &
            units="m s-1", long_name="Basal velocity in x", &
            standard_name="land_ice_basal_x_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "yvelbase", ylmo%dyn%now%uy_b / yr_to_sec, &
            units="m s-1", long_name="Basal velocity in y", &
            standard_name="land_ice_basal_y_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "zvelbase", uz_b_masked, &
            units="m s-1", long_name="Basal velocity in z", &
            standard_name="land_ice_basal_upward_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "xvelmean", ux_aa / yr_to_sec, &
            units="m s-1", long_name="Mean velocity in x", &
            standard_name="land_ice_vertical_mean_x_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "yvelmean", uy_aa / yr_to_sec, &
            units="m s-1", long_name="Mean velocity in y", &
            standard_name="land_ice_vertical_mean_y_velocity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litemptop", T_top_ice, &
            units="K", long_name="Surface temperature", &
            standard_name="temperature_at_top_of_ice_sheet_model", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litempavg", T_avg, &
            units="K", long_name="Depth-averaged ice temperature", &
            standard_name="land_ice_temperature", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litempbotgr", T_base_grnd, &
            units="K", long_name="Basal temperature beneath grounded ice sheet", &
            standard_name="temperature_at_base_of_ice_sheet_model", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litempbotfl", T_base_flt, &
            units="K", long_name="Basal temperature beneath floating ice shelf", &
            standard_name="temperature_at_base_of_ice_sheet_model", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litempgradgr", dTdz_base_grnd, &
            units="K m-1", long_name="Vertical basal temperature gradient beneath grounded ice sheet", &
            standard_name="", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "litempgradfl", dTdz_base_flt, &
            units="K m-1", long_name="Vertical basal temperature gradient beneath floating ice shelf", &
            standard_name="", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "strbasemag", ylmo%dyn%now%taub, &
            units="Pa", long_name="Basal drag", &
            standard_name="land_ice_basal_drag", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "sftgif", ylmo%tpo%now%f_ice, &
            units="1", long_name="Land ice area fraction", &
            standard_name="land_ice_area_fraction", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "sftgrf", ylmo%tpo%now%f_grnd, &
            units="1", long_name="Grounded ice sheet area fraction", &
            standard_name="grounded_ice_sheet_area_fraction", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "sftflf", MAX(ylmo%tpo%now%f_ice - ylmo%tpo%now%f_grnd, 0.0_wp), &
            units="1", long_name="Floating ice sheet area fraction", &
            standard_name="floating_ice_shelf_area_fraction", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        ! ====================================================================
        ! 2-D FL variables  (yearly-average flux)
        ! ====================================================================
        
        ! Time-invariant in most setups; written without time dimension.
        call nc_write(filename, "hfgeoubed", ylmo%bnd%Q_geo * 1.0e3_wp, &
            units="W m-2", long_name="Geothermal heat flux", &
            standard_name="upward_geothermal_heat_flux_in_land_ice", &
            dim1="xc", dim2="yc", start=[1,1], ncid=ncid)
        
        call nc_write(filename, "acabf", ylmo%tpo%now%smb * esm_correction, &
            units="kg m-2 s-1", long_name="Surface mass balance flux", &
            standard_name="land_ice_surface_specific_mass_balance_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "libmassbfgr", bmb_grnd_masked * esm_correction, &
            units="kg m-2 s-1", long_name="Basal mass balance flux beneath grounded ice", &
            standard_name="land_ice_basal_specific_mass_balance_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "libmassbffl", bmb_shlf_masked * esm_correction, &
            units="kg m-2 s-1", long_name="Basal mass balance flux beneath floating ice", &
            standard_name="land_ice_basal_specific_mass_balance_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        call nc_write(filename, "dlithkdt", ylmo%tpo%now%dHidt / yr_to_sec, &
            units="m s-1", long_name="Ice thickness imbalance", &
            standard_name="tendency_of_land_ice_thickness", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        ! Kinematic flux through ice-front cells (mask_frnt==1, floating).
        ! Uses same formula as CalvingMIP
        call nc_write(filename, "licalvf", flux_clv_2d, &
            units="kg m-2 s-1", long_name="Calving flux", &
            standard_name="land_ice_specific_mass_flux_due_to_calving", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        ! ligroundf : Grounding-line flux                          [MANDATORY]
        call nc_write(filename, "ligroundf", flux_grl_2d, &
            units="kg m-2 s-1", long_name="Grounding line flux", &
            standard_name="land_ice_specific_grounding_line_flux", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        ! tfbase : Thermal forcing at ice base, floating           [optional]
        call nc_write(filename, "tfbase", tfbase, &
            units="K", long_name="Thermal forcing at the ice base", &
            standard_name="", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)
        
        ! ---- close ----------------------------------------------------------
        call nc_close(ncid)
        
        return
        
    end subroutine write_step_2D_cmip
        
        
    subroutine write_step_1D_cmip(ylmo, mshlf, filename, time)
        ! Writes all mandatory scalar ISMIP7 (ISM_2026) variables.
        ! Flux computation follows yelmo_calving.f90: kinematic uxy_bar*H_ice*rho_ice.
        ! -------------------------------------------------------------------------
        
        implicit none
        
        type(yelmo_class),    intent(IN) :: ylmo
        type(marshelf_class), intent(IN) :: mshlf
        character(len=*),     intent(IN) :: filename
        real(wp),             intent(IN) :: time
        
        ! ---- local variables ------------------------------------------------
        type(yregions_class) :: reg
        
        integer  :: ncid, n
        real(wp) :: rho_ice, esm_correction, yr_to_sec
        real(wp) :: dx, dy
        
        real(wp) :: smb_tot          ! total SMB           [m3 yr-1]
        real(wp) :: bmb_grnd_tot     ! total BMB grounded  [m3 yr-1]
        real(wp) :: bmb_shlf_t       ! total BMB floating  [m3 yr-1]
        real(wp) :: flux_grl         ! total GL flux       [kg yr-1]
        real(wp) :: flux_clv         ! total calving flux  [kg yr-1]
        
        logical, allocatable :: mask_tot(:,:)
        logical, allocatable :: mask_grnd(:,:)
        logical, allocatable :: mask_flt(:,:)
        logical, allocatable :: mask_grl(:,:)    ! grounding-line cells
        logical, allocatable :: mask_frnt(:,:)   ! ice-front cells
        
        ! ---- allocate -------------------------------------------------------
        allocate(mask_tot  (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_grnd (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_flt  (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_grl  (ylmo%grd%G%nx, ylmo%grd%G%ny))
        allocate(mask_frnt (ylmo%grd%G%nx, ylmo%grd%G%ny))
        
        ! === Unit conversion factors =========================================
        ! Taken from Yelmo's own constants, so the output cannot disagree with
        ! the model that produced it. esm_correction converts a rate in
        ! [<unit> yr-1] of ice to [kg ... s-1]; it collapses to rho_ice/sec_year,
        ! since the old two-step form was (rho_w/sec_year)*(rho_ice/rho_w) with
        ! rho_w/sec_year rounded to 3.2e-5 (1 % high) and sec_year written as
        ! 31556952 (the Gregorian year) rather than Yelmo's own.
        rho_ice        = ylmo%bnd%c%rho_ice
        yr_to_sec      = ylmo%bnd%c%sec_year
        esm_correction = rho_ice / yr_to_sec
        
        dx = ylmo%grd%G%dx
        dy = ylmo%grd%G%dy
        
        ! ---- masks (updated to match yelmo_calving.f90) ---------------------
        mask_tot  = (ylmo%tpo%now%H_ice .gt. 0.0_wp)
        mask_grnd = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .gt. 0.0_wp)
        mask_flt  = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp)
        
        ! Grounding-line cells: grounded ice where mask_grz == 0
        mask_grl  = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .gt. 0.0_wp &
                        .and. ylmo%tpo%now%mask_grz .eq. 0.0_wp)
        
        ! Ice-front cells: floating ice where mask_frnt == 1
        mask_frnt = (ylmo%tpo%now%H_ice .gt. 0.0_wp .and. ylmo%tpo%now%f_grnd .eq. 0.0_wp &
                        .and. ylmo%tpo%now%mask_frnt .eq. 1.0_wp)
        
        ! ---- regional object (pre-computed by Yelmo) ------------------------
        reg = ylmo%reg
        
        ! ---- integrated fluxes ----------------------------------------------
        ! SMB and BMB: area-integrated [m3 yr-1], converted via esm_correction
        smb_tot      = sum(ylmo%bnd%smb,       mask=mask_tot)  * (dx * dy)
        bmb_grnd_tot = sum(ylmo%tpo%now%bmb,   mask=mask_grnd) * (dx * dy)
        bmb_shlf_t   = sum(ylmo%tpo%now%bmb,   mask=mask_flt)  * (dx * dy)
        
        ! Grounding-line flux: kinematic, [kg yr-1]
        ! Matches yelmo_calving.f90: uxy_bar * H_ice * rho_ice * dx
        if (count(mask_grl) .gt. 0) then
            flux_grl = sum(ylmo%dyn%now%uxy_bar * ylmo%tpo%now%H_ice * rho_ice, mask=mask_grl) * dx
        else
            flux_grl = 0.0_wp
        end if
        
        ! Calving flux: sum of cmb_flt and cmb_grnd, [kg yr-1]
        flux_clv = sum(ylmo%tpo%now%cmb_flt + ylmo%tpo%now%cmb_grnd) * (dx * dy)  ! [m3 yr-1]
        
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
        call nc_write(filename, "tendacabf", smb_tot * esm_correction, &
            units="kg s-1", long_name="Total SMB flux", &
            standard_name="tendency_of_land_ice_mass_due_to_surface_mass_balance", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendlibmassbfgr : Total BMB flux, grounded               [MANDATORY]
        call nc_write(filename, "tendlibmassbfgr", bmb_grnd_tot * esm_correction, &
            units="kg s-1", long_name="Total BMB flux beneath grounded ice", &
            standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendlibmassbffl : Total BMB flux, floating               [MANDATORY]
        call nc_write(filename, "tendlibmassbffl", bmb_shlf_t * esm_correction, &
            units="kg s-1", long_name="Total BMB flux beneath floating ice", &
            standard_name="tendency_of_land_ice_mass_due_to_basal_mass_balance", &
            dim1="time", start=[n], ncid=ncid)
        
        ! tendlicalvf : Total calving flux                         [MANDATORY]
        call nc_write(filename, "tendlicalvf", flux_clv * esm_correction, &
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
        ! Kinematic, convert from kg yr-1 to kg s-1.
        call nc_write(filename, "tendligroundf", flux_grl / yr_to_sec, &
            units="kg s-1", long_name="Total grounding line flux", &
            standard_name="tendency_of_grounded_ice_sheet_mass", &
            dim1="time", start=[n], ncid=ncid)
        
        ! ---- close ----------------------------------------------------------
        call nc_close(ncid)
        
        return
        
    end subroutine write_step_1D_cmip

end module cmip_output
