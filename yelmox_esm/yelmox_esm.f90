program yelmox_esm
    ! Multigrid yelmox driver with ESM climatic forcing (single domain).
    !
    ! Reuses the shared multigrid domain machinery (libs/yelmox_domain.f90) for
    ! everything except climate: one ice_domain with each sub-model on its own
    ! configurable grid, the hi-res topography hub, and the coupler maps. The
    ! climate/ocean forcing comes from the ESM module (libs/esm.f90) rather than
    ! snapclim, so the driver owns an esm_forcing_class and calls its own
    ! step_climate_esm / step_marine_shelf_esm in place of the snapclim-based
    ! step_climate / step_marine_shelf.
    !
    ! esm is a first-class multigrid component: it runs entirely on its own grid
    ! (grid_clim, the climate grid) and each output is remapped to the consumer
    ! module's grid at coupling time, exactly like snapclim -- atmosphere fields to
    ! grid_smb (smbpal), the depth-interpolated ocean forcing to grid_mshlf
    ! (marine_shelf). Geometry comes from the hi-res hub, remapped to whichever grid
    ! a step needs. marshelf_interp_shelf reads only mshlf%par (grid-agnostic), so
    ! the ESM ocean interpolation runs on grid_clim and the resulting T_shlf/S_shlf
    ! are remapped to grid_mshlf. With grid_clim == grid_smb == grid_mshlf ==
    ! grid_yelmo == grid_name every remap is an identity copy, reproducing
    ! yelmox_esm.f90; set grid_clim to a coarse ESM grid and it genuinely fans out.
    !
    ! Output (yelmo / yelmo_sm / yelmo_ts_esm and the CMIP-formatted files) is kept
    ! identical to yelmox_esm.f90 via yelmox_esm_output. See docs/multigrid.md.

    use nml
    use ncio
    use timestepping
    use timeout
    use yelmo
    use ice_optimization

    use esm_forcing
    use fastisostasy    ! also reexports barysealevel (bsl_*)
    use marine_shelf, only : marshelf_update, marshelf_interp_shelf, ocn_variable_extrapolation
    use smbpal,       only : smbpal_update_monthly, smbpal_update_monthly_equil

    use yelmox_domain
    use yelmox_esm_output

    implicit none

    character(len=512) :: path_par, esm_path_par
    character(len=512) :: outfldr
    character(len=512) :: file2D, file2Dsm, file1D_esm, file1D_cmip, file2D_cmip
    character(len=512) :: file_isos, file_bsl

    type(tstep_class)   :: ts
    type(ice_domain)    :: dom
    type(bsl_class)     :: bsl          ! shared, driver-owned barystatic sea level
    type(esm_forcing_class) :: esm      ! driver-owned climate forcing (replaces snapclim)
    type(timeout_class) :: tm_1D, tm_2D, tm_2Dsm

    ! ESM run control (esm-specific config; the shared domain_ctl carries the
    ! rest, and the timeline is driver-owned via timeline_init)
    type esm_ctl_params
        character(len=56) :: run_step
        real(wp) :: dtt, time_equil
        real(wp) :: time_ref(2), time_hist(2), time_proj(2), time_esm_ref(2)
        character(len=56) :: clim_var
        integer  :: clim_seed
        logical  :: kill_shelves
        logical  :: use_esm, use_smb, use_var, use_proj, use_hist
        logical  :: write_formatted
        real(wp) :: dt_formatted
        character(len=512) :: par_file
        character(len=56)  :: experiment, esm_name
    end type
    type(esm_ctl_params) :: ec

    integer :: seed_size
    integer, allocatable :: seed(:)

    ! ================= SETUP ==================================================

    call yelmo_load_command_line_args(path_par)

    ! ESM run control: [ctrl] run_step selects the run phase, [esm] holds the
    ! experiment identity + physics parameters, and the [run_step] group the
    ! phase's esm timing (the timeline keys of the same group are read by
    ! timeline_init / domain_init below).
    call esm_ctl_load(ec, esm, path_par)

    ! Seed the RNG for climate-variability reproducibility.
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = ec%clim_seed
    call random_seed(put=seed)

    ! Single-domain runs write to the run dir.
    outfldr     = "./"
    ! File naming convention: 2D output -> "yelmo", 1D timeseries -> "yelmo_ts".
    file2D      = trim(outfldr)//"yelmo.nc"
    file2Dsm    = trim(outfldr)//"yelmo_sm.nc"
    file1D_esm  = trim(outfldr)//"yelmo_ts_esm.nc"
    file1D_cmip = trim(outfldr)//"yelmo_ts_cmip.nc"
    file2D_cmip = trim(outfldr)//"yelmo_cmip.nc"
    file_isos   = trim(outfldr)//"fastisostasy.nc"
    file_bsl    = trim(outfldr)//"bsl.nc"

    ! Timestepping (driver-owned; the [run_step] group holds this run phase's
    ! timeline, with tstep_const applied as a calendar constant).
    call timeline_init(ts, ec%dtt, path_par, trim(ec%run_step), &
                       time_ref=2000.0_wp, cal=.true.)

    write(*,*)
    write(*,*) "yelmox_esm: run_step = "//trim(ec%run_step)
    write(*,*) "  time_init/end/dtt: ", ts%time_init, ts%time_end, ec%dtt
    write(*,*) "  esm experiment:    "//trim(ec%experiment)//"  use_esm = ", ec%use_esm
    write(*,*)

    ! Shared, driver-owned barystatic sea level.
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, ts%time_rel)

    ! Initialize the domain (sub-models + hub + coupler maps), skipping snapclim:
    ! this driver supplies its own climate via the esm object. The domain reads
    ! the timeline values it needs from the [run_step] group.
    call domain_init(dom, path_par, ts%time, init_climate=.false., &
                     timeline_group=trim(ec%run_step))

    ! Regions of interest for 1D output (must precede the first yelmo_update).
    call domain_regions_init(dom, trim(outfldr))

    ! Initialize the ESM forcing on its own grid (grid_clim, the climate grid).
    ! esm runs entirely on this grid; its outputs are remapped to the consumer
    ! grids (grid_smb for smbpal, grid_mshlf for marine_shelf) at coupling time.
    esm_path_par = trim(outfldr)//"/"//trim(ec%par_file)
    call esm_forcing_init(esm, esm_path_par, dom%ctl%domain, trim(dom%ctl%grid_clim), &
                          run_type=ec%run_step, gcm=ec%esm_name, experiment=ec%experiment, &
                          use_esm=ec%use_esm, use_smb=ec%use_smb, use_var=ec%use_var, &
                          use_hist=ec%use_hist, use_proj=ec%use_proj)

    write(*,*)
    write(*,*) "yelmox_esm: domain initialized"
    write(*,*) "  domain      : "//trim(dom%ctl%domain)
    write(*,*) "  Yelmo grid  : "//trim(dom%ctl%grid_yelmo), dom%yelmo%grd%G%nx, dom%yelmo%grd%G%ny
    write(*,*) "  esm  grid   : "//trim(dom%ctl%grid_clim)
    write(*,*) "  smb  grid   : "//trim(dom%ctl%grid_smb)
    write(*,*) "  mshlf grid  : "//trim(dom%ctl%grid_mshlf)
    write(*,*) "  coupler maps: ", dom%cpl%nmaps
    write(*,*)

    ! Sea level + isostasy output init (parity with yelmox_esm.f90; init only).
    call bsl_write_init(bsl, file_bsl, ts%time)
    call isos_write_init_extended(dom%isos, file_isos, ts%time)

    ! Cold start: build the initial boundary state. Restart: restore the bundle
    ! (incl. the shared bsl) and rebuild the hub, then re-establish the forcing.
    if (trim(dom%ctl%restart) == "None") then
        call esm_cold_start(dom, esm, ec, ts, bsl)
    else
        call domain_startup(dom, ts, bsl)
        call step_climate_esm(dom, esm, ec, ts)
        call step_marine_shelf_esm(dom, esm, ec, ts)
    end if

    ! ================= OUTPUT SETUP ==========================================

    call timeout_init(tm_1D,   path_par, "tm_1D",   "small",  ts%time_init, ts%time_end)
    call timeout_init(tm_2D,   path_par, "tm_2D",   "heavy",  ts%time_init, ts%time_end)
    call timeout_init(tm_2Dsm, path_par, "tm_2Dsm", "medium", ts%time_init, ts%time_end)

    call yelmo_write_init(dom%yelmo, file2D,   time_init=ts%time, units="years")
    call yelmo_write_init(dom%yelmo, file2Dsm, time_init=ts%time, units="years")
    call yelmo_regions_write(dom%yelmo, ts%time, init=.true., units="years")
    call yelmo_write_reg_init(dom%yelmo, file1D_esm, time_init=ts%time, units="years", &
                              mask=(dom%yelmo%bnd%mask_ice /= MASK_ICE_NONE))

    if (ec%write_formatted) then
        call yelmo_write_init(dom%yelmo, file2D_cmip, time_init=ts%time, units="years")
        call yelmo_write_reg_init(dom%yelmo, file1D_cmip, time_init=ts%time, units="years", &
                                  mask=(dom%yelmo%bnd%mask_ice /= MASK_ICE_NONE))
    end if

    ! ================= MAIN TIME LOOP ========================================

    call tstep_print_header(ts)
    do while (.not. ts%is_finished)

        call tstep_update(ts, dom%ctl%dtt)
        call tstep_print(ts)

        ! Shared sea level: once per step, before the domain advances.
        call bsl_update(bsl, ts%time_rel)

        ! Shared multigrid primitives (isostasy / ice sheet / hub / optimization).
        ! The ESM-owned Yelmo input (Qd) is landed just before the ice sheet runs,
        ! alongside the shared couplers invoked inside step_icesheet.
        call step_optimize(dom, ts)
        call step_isostasy(dom, ts, bsl)
        call couple_esm_extras_to_yelmo(dom, esm)
        call step_icesheet(dom, ts)
        call refresh_htopo(dom)

        ! ESM climate + ocean forcing (driver-owned, replaces snapclim steps).
        call step_climate_esm(dom, esm, ec, ts)
        call step_marine_shelf_esm(dom, esm, ec, ts)

        ! === Model output ===
        if (timeout_check(tm_2Dsm, ts%time)) &
            call write_step_2D_small(dom%yelmo, dom%isos, esm, dom%mshlf, dom%smb, &
                                     ec%use_smb, file2Dsm, ts%time)
        if (timeout_check(tm_2D, ts%time)) &
            call write_step_2D_combined(dom%yelmo, dom%isos, esm, dom%mshlf, dom%smb, &
                                        ec%use_smb, file2D, ts%time)
        if (timeout_check(tm_1D, ts%time)) then
            call yelmo_regions_write(dom%yelmo, ts%time)
            call write_1D_esm(dom%yelmo, esm, dom%mshlf, file1D_esm, ts%time)
        end if

        if (ec%write_formatted) then
            if (tstep_due(ts%time_elapsed, ec%dt_formatted)) then
                call write_step_2D_cmip(dom%yelmo, dom%mshlf, file2D_cmip, ts%time)
                call write_step_1D_cmip(dom%yelmo, dom%mshlf, file1D_cmip, ts%time)
            end if
        end if

        ! === Restart bundle (domain + shared bsl) ===
        if (tstep_due(ts%time, dom%ctl%dt_restart)) then
            call run_restart_write(dom, bsl, ts%time)
        end if

    end do

    ! Final state + restart bundle (incl. shared bsl).
    call write_step_2D_combined(dom%yelmo, dom%isos, esm, dom%mshlf, dom%smb, &
                                ec%use_smb, file2D, ts%time)
    call yelmo_regions_write(dom%yelmo, ts%time)
    call run_restart_write(dom, bsl, ts%time)

    write(*,*)
    write(*,*) "yelmox_esm: run complete at time =", ts%time
    write(*,*) "  H_ice max =", maxval(dom%yelmo%tpo%now%H_ice)

    call yelmo_end(dom%yelmo, time=ts%time)
    deallocate(seed)

contains

    ! ---------------------------------------------------------------------------
    subroutine esm_ctl_load(ec, esm, path_par)
        ! Load the ESM run control: [ctrl] run_step, the [esm] experiment
        ! identity + physics parameters (the latter straight into the esm
        ! object), and the selected [run_step] group's esm timing periods.
        type(esm_ctl_params),    intent(inout) :: ec
        type(esm_forcing_class), intent(inout) :: esm
        character(len=*),        intent(in)    :: path_par

        call nml_read(path_par, "ctrl", "run_step", ec%run_step)

        ! [esm] group: experiment identity + physics parameters.
        call nml_read(path_par, "esm", "par_file",        ec%par_file)
        call nml_read(path_par, "esm", "experiment",      ec%experiment)
        call nml_read(path_par, "esm", "esm_name",        ec%esm_name)
        call nml_read(path_par, "esm", "use_esm",         ec%use_esm)
        call nml_read(path_par, "esm", "use_smb",         ec%use_smb)
        call nml_read(path_par, "esm", "use_var",         ec%use_var)
        call nml_read(path_par, "esm", "use_proj",        ec%use_proj)
        call nml_read(path_par, "esm", "use_hist",        ec%use_hist)
        call nml_read(path_par, "esm", "write_formatted", ec%write_formatted)
        call nml_read(path_par, "esm", "dt_formatted",    ec%dt_formatted)
        call nml_read(path_par, "esm", "lapse",           esm%lapse)
        call nml_read(path_par, "esm", "f_p",             esm%beta_p)
        call nml_read(path_par, "esm", "f_ocn",           esm%f_ocn)
        call nml_read(path_par, "esm", "f_polar",         esm%f_polar)
        call nml_read(path_par, "esm", "dT_threshold",    esm%dT_lim)
        call nml_read(path_par, "esm", "grid_src",        esm%grid_src)

        ! [run_step] group: esm reference/history/projection periods + switches.
        call nml_read(path_par, trim(ec%run_step), "time_equil",   ec%time_equil)
        call nml_read(path_par, trim(ec%run_step), "time_ref",     ec%time_ref)
        call nml_read(path_par, trim(ec%run_step), "time_hist",    ec%time_hist)
        call nml_read(path_par, trim(ec%run_step), "time_proj",    ec%time_proj)
        call nml_read(path_par, trim(ec%run_step), "time_esm_ref", ec%time_esm_ref)
        call nml_read(path_par, trim(ec%run_step), "clim_var",     ec%clim_var)
        call nml_read(path_par, trim(ec%run_step), "clim_seed",    ec%clim_seed)
        call nml_read(path_par, trim(ec%run_step), "kill_shelves", ec%kill_shelves)
    end subroutine esm_ctl_load

    ! ---------------------------------------------------------------------------
    subroutine esm_cold_start(dom, esm, ec, ts, bsl)
        ! Build the initial boundary state for a cold start (no restart), mirroring
        ! yelmox_esm.f90's initialization: isostasy reference/state, first climate
        ! forcing, snowpack equilibration, Yelmo state init, optional shelf kill,
        ! and (for the "opt" spinup) the Yelmo equilibration passes.
        type(ice_domain),     intent(inout) :: dom
        type(esm_forcing_class), intent(inout) :: esm
        type(esm_ctl_params), intent(in)    :: ec
        type(tstep_class),    intent(in)    :: ts
        type(bsl_class),      intent(inout) :: bsl

        real(wp), allocatable :: z_bed_ref_i(:,:), H_ice_ref_i(:,:)
        real(wp), allocatable :: z_bed_i(:,:), H_ice_i(:,:)
        character(len=256) :: gi, gy

        gi = trim(dom%ctl%grid_isos)
        gy = trim(dom%ctl%grid_yelmo)

        ! Isostasy reference + initial state (isostasy runs on grid_isos).
        call remap(dom, dom%yelmo%bnd%z_bed_ref, gy, z_bed_ref_i, gi, "bilin")
        call remap(dom, dom%yelmo%bnd%H_ice_ref, gy, H_ice_ref_i, gi, "bilin")
        call isos_init_ref(dom%isos, z_bed_ref_i, H_ice_ref_i)
        call remap(dom, dom%yelmo%bnd%z_bed,     gy, z_bed_i, gi, "bilin")
        call remap(dom, dom%yelmo%tpo%now%H_ice, gy, H_ice_i, gi, "bilin")
        call isos_init_state(dom%isos, z_bed_i, H_ice_i, ts%time, bsl)
        call couple_isostasy_to_yelmo(dom)

        ! Hub mirror + first ESM forcing (init=.true. runs the smbpal ITM equil).
        call refresh_htopo(dom)
        call step_climate_esm(dom, esm, ec, ts, init=.true.)
        call step_marine_shelf_esm(dom, esm, ec, ts)

        ! Assemble the Yelmo boundary state from the freshly produced outputs.
        call couple_smb_to_yelmo(dom)
        call couple_marine_to_yelmo(dom)
        call couple_esm_extras_to_yelmo(dom, esm)

        ! Cold-start friction guess for the optimization.
        if (trim(dom%ctl%equil_method) == "opt") &
            dom%yelmo%dyn%now%cb_ref = dom%opt%cf_init

        ! Initialize Yelmo state variables (cold base).
        call yelmo_print_bound(dom%yelmo%bnd)
        call yelmo_init_state(dom%yelmo, time=ts%time, thrm_method="robin-cold")

        ! Optional: kill ice shelves beyond present-day extent.
        if (ec%kill_shelves) then
            where(dom%yelmo%dta%pd%mask_bed .eq. mask_bed_ocean) &
                dom%yelmo%bnd%mask_ice = MASK_ICE_NONE
        end if

        ! Spinup ("opt") equilibration passes (cold start only).
        if (trim(ec%run_step) == "spinup" .and. dom%ctl%with_ice_sheet) then
            call yelmo_update_equil(dom%yelmo, ts%time, time_tot=1.0_wp, dt=1.0_wp, topo_fixed=.false.)
            if (trim(dom%ctl%equil_method) == "opt" .and. ec%time_equil > 0.0_wp) &
                call yelmo_update_equil(dom%yelmo, ts%time, time_tot=ec%time_equil, &
                                        dt=ec%dtt, topo_fixed=.true.)
        end if

    end subroutine esm_cold_start

    ! ---------------------------------------------------------------------------
    subroutine step_climate_esm(dom, esm, ec, ts, init)
        ! ESM climate + surface mass balance. esm runs entirely on its own grid
        ! (grid_clim): geometry is remapped from the hub onto it, the three esm
        ! updates produce the atmospheric/ocean anomaly fields there, and the
        ! atmospheric forcing is remapped to grid_smb for smbpal (or the direct-SMB
        ! path), with smb / T_srf aggregated back to the Yelmo grid. The ocean
        ! anomalies (esm%dto/dso, esm%to_ref/so_ref) stay on grid_clim for
        ! step_marine_shelf_esm.
        type(ice_domain),        intent(inout) :: dom
        type(esm_forcing_class), intent(inout) :: esm
        type(esm_ctl_params),    intent(in)    :: ec
        type(tstep_class),       intent(in)    :: ts
        logical, intent(in), optional          :: init

        real(wp), allocatable :: z_srf_e(:,:), H_ice_e(:,:), z_bed_e(:,:)
        real(wp), allocatable :: f_grnd_e(:,:), z_sl_e(:,:), basins_e(:,:)
        real(wp), allocatable :: tas_s(:,:,:), pr_s(:,:,:), z_srf_s(:,:), H_ice_s(:,:)
        real(wp), allocatable :: smb_ann_s(:,:), dsmb_s(:,:), dsmbdz_s(:,:)
        real(wp), allocatable :: pd_zsrf_s(:,:), tsrf_s(:,:)
        character(len=256) :: ge, gn, gy, gs
        logical :: is_init

        is_init = .false.
        if (present(init)) is_init = init

        ge = trim(dom%ctl%grid_clim)    ! esm grid (esm's own working grid)
        gn = trim(dom%ctl%grid_name)    ! hub grid
        gy = trim(dom%ctl%grid_yelmo)
        gs = trim(dom%ctl%grid_smb)

        ! Geometry: hub -> esm grid.
        call remap(dom, dom%topo%z_srf,   gn, z_srf_e,  ge, "bilin")
        call remap(dom, dom%topo%H_ice,   gn, H_ice_e,  ge, "bilin")
        call remap(dom, dom%topo%z_bed,   gn, z_bed_e,  ge, "bilin")
        call remap(dom, dom%topo%f_grnd,  gn, f_grnd_e, ge, "bilin")
        call remap(dom, dom%topo%z_sl,    gn, z_sl_e,   ge, "bilin")
        call remap(dom, dom%topo%basins,  gn, basins_e, ge, "nn")

        ! Step 1: reference climatology (lapse-rate / precip scaling to z_srf).
        call esm_clim_update(esm, z_srf_e, ts%time, ec%time_ref, ec%use_smb, &
                             dom%ctl%domain, ge)

        ! Extrapolate reference ocean into ice-shelf interiors.
        if (dom%mshlf%par%extrap_shlf) then
            call ocn_variable_extrapolation(esm%to_ref%var(:,:,:,1), H_ice_e, basins_e, &
                                            -esm%to_var_ref%z, z_bed_e)
            call ocn_variable_extrapolation(esm%so_ref%var(:,:,:,1), H_ice_e, basins_e, &
                                            -esm%so_var_ref%z, z_bed_e)
        end if

        ! Step 2: anomaly fields (historical / projection / homogeneous). mshlf is
        ! passed for its params only (marshelf_interp_shelf is grid-agnostic), so
        ! the ocean anomalies dto/dso are produced on the esm grid.
        call esm_forcing_update(esm, dom%mshlf, ts%time, ec%use_esm, ec%time_ref, &
                                ec%time_hist, ec%time_proj, ec%time_esm_ref, dom%ctl%domain, &
                                H_ice_e, basins_e, z_bed_e, f_grnd_e, z_sl_e, ec%use_smb, &
                                use_ref_atm=.false., use_ref_ocn=.false.)

        ! Step 3: variability anomaly.
        call esm_variability_update(esm, dom%mshlf, ts%time, ec%dtt, ec%clim_var, ec%time_ref, &
                                    H_ice_e, basins_e, z_bed_e, f_grnd_e, z_sl_e, ec%use_var, &
                                    use_ref_atm=.false., use_ref_ocn=.false.)

        ! === Atmospheric boundary conditions (SMB on grid_smb) ===
        if (ec%use_smb) then
            ! Direct SMB from the esm fields: remap each term esm grid -> grid_smb,
            ! evaluate the elevation correction there against the (remapped) surface.
            call remap(dom, esm%smb_ann,               ge, smb_ann_s, gs, "bilin")
            call remap(dom, sum(esm%dsmb, dim=3)/12.0_wp, ge, dsmb_s,  gs, "bilin")
            call remap(dom, esm%dsmbdz,                ge, dsmbdz_s,  gs, "bilin")
            call remap(dom, dom%yelmo%dta%pd%z_srf,    gy, pd_zsrf_s, gs, "bilin")
            call remap(dom, dom%topo%z_srf,            gn, z_srf_s,   gs, "bilin")
            call remap(dom, dom%topo%H_ice,            gn, H_ice_s,   gs, "bilin")
            call remap(dom, sum(esm%t2m + esm%dts + esm%dts_var, dim=3)/12.0_wp, &
                       ge, tsrf_s, gs, "bilin")
            dom%smb%ann%smb  = smb_ann_s + dsmb_s - dsmbdz_s*(pd_zsrf_s - z_srf_s)
            dom%smb%ann%tsrf = tsrf_s
            where(H_ice_s > 0.0_wp .and. tsrf_s > 273.15_wp) dom%smb%ann%tsrf = 273.15_wp
        else
            ! smbpal on grid_smb: atmospheric forcing from the esm grid, geometry
            ! from the hub. init=.true. runs the ITM snowpack equilibration first.
            call remap(dom, esm%t2m + esm%dts + esm%dts_var, ge, tas_s, gs, "bilin")
            call remap(dom, esm%pr  * esm%dpr * esm%dpr_var, ge, pr_s,  gs, "bilin")
            call remap(dom, dom%topo%z_srf, gn, z_srf_s, gs, "bilin")
            call remap(dom, dom%topo%H_ice, gn, H_ice_s, gs, "bilin")
            if (is_init .and. trim(dom%smb%par%abl_method) == "itm") then
                call remap(dom, esm%t2m + esm%dts, ge, tas_s, gs, "bilin")
                call remap(dom, esm%pr  * esm%dpr, ge, pr_s,  gs, "bilin")
                call smbpal_update_monthly_equil(dom%smb, tas_s, pr_s, z_srf_s, H_ice_s, &
                                                 ts%time_rel, time_equil=100.0_wp)
                call remap(dom, esm%t2m + esm%dts + esm%dts_var, ge, tas_s, gs, "bilin")
                call remap(dom, esm%pr  * esm%dpr * esm%dpr_var, ge, pr_s,  gs, "bilin")
            end if
            call smbpal_update_monthly(dom%smb, tas_s, pr_s, z_srf_s, H_ice_s, ts%time_rel)
        end if

        ! smb / tsrf now live on grid_smb in dom%smb%ann; couple_smb_to_yelmo lands
        ! them on the Yelmo grid. Subglacial discharge (esm%Qd_ann) is landed by
        ! couple_esm_extras_to_yelmo. Both run in the ice-sheet coupling step.

    end subroutine step_climate_esm

    ! ---------------------------------------------------------------------------
    subroutine step_marine_shelf_esm(dom, esm, ec, ts)
        ! ESM ocean boundary conditions + marine-shelf basal melt. The ocean forcing
        ! is produced on the esm grid (grid_clim): interpolate the reference ocean to
        ! shelf depth (marshelf_interp_shelf reads only mshlf%par, so it is
        ! grid-agnostic) and add the esm ocean anomalies from step_climate_esm. The
        ! resulting T_shlf/S_shlf are remapped to grid_mshlf, where marshelf_update
        ! runs against the hub geometry, and bmb_shlf/T_shlf are aggregated back to
        ! the Yelmo grid.
        type(ice_domain),        intent(inout) :: dom
        type(esm_forcing_class), intent(inout) :: esm
        type(esm_ctl_params),    intent(in)    :: ec
        type(tstep_class),       intent(in)    :: ts

        real(wp), allocatable :: H_ice_e(:,:), z_bed_e(:,:), f_grnd_e(:,:), z_sl_e(:,:)
        real(wp), allocatable :: T_shlf_e(:,:), S_shlf_e(:,:), dT_shlf_e(:,:)
        real(wp), allocatable :: H_ice_m(:,:), z_bed_m(:,:), f_grnd_m(:,:), z_sl_m(:,:)
        real(wp), allocatable :: regions_m(:,:), basins_m(:,:)
        character(len=256) :: ge, gm, gn

        if (.not. dom%ctl%with_marine_shelf) return

        ge = trim(dom%ctl%grid_clim)    ! esm grid (ocean fields live here)
        gm = trim(dom%ctl%grid_mshlf)   ! marine-shelf grid
        gn = trim(dom%ctl%grid_name)    ! hub grid

        ! --- Ocean forcing on the esm grid ---
        ! Geometry for the shelf-depth interpolation: hub -> esm grid.
        call remap(dom, dom%topo%H_ice,  gn, H_ice_e,  ge, "bilin")
        call remap(dom, dom%topo%z_bed,  gn, z_bed_e,  ge, "bilin")
        call remap(dom, dom%topo%f_grnd, gn, f_grnd_e, ge, "bilin")
        call remap(dom, dom%topo%z_sl,   gn, z_sl_e,   ge, "bilin")

        allocate(T_shlf_e(size(H_ice_e,1), size(H_ice_e,2)))
        allocate(S_shlf_e(size(H_ice_e,1), size(H_ice_e,2)))
        call marshelf_interp_shelf(T_shlf_e, dom%mshlf, esm%to_ref%var(:,:,:,1), &
                                   H_ice_e, z_bed_e, f_grnd_e, z_sl_e, -esm%to_ref%z)
        call marshelf_interp_shelf(S_shlf_e, dom%mshlf, esm%so_ref%var(:,:,:,1), &
                                   H_ice_e, z_bed_e, f_grnd_e, z_sl_e, -esm%so_ref%z)
        T_shlf_e = T_shlf_e + esm%dto + esm%dto_var
        S_shlf_e = S_shlf_e + esm%dso + esm%dso_var

        ! Send the ocean forcing to the marine-shelf grid (esm grid -> grid_mshlf).
        call remap(dom, T_shlf_e, ge, dom%mshlf%now%T_shlf, gm, "bilin")
        call remap(dom, S_shlf_e, ge, dom%mshlf%now%S_shlf, gm, "bilin")

        if (trim(dom%ctl%domain) == "Greenland") then
            dT_shlf_e = T_shlf_e + esm%dto
            call remap(dom, dT_shlf_e, ge, dom%mshlf%now%dT_shlf, gm, "bilin")
            dom%mshlf%par%tf_method = 2
        end if

        ! --- Marine-shelf basal melt on grid_mshlf ---
        ! Geometry + masks: hub -> mshlf grid.
        call remap(dom, dom%topo%H_ice,   gn, H_ice_m,   gm, "bilin")
        call remap(dom, dom%topo%z_bed,   gn, z_bed_m,   gm, "bilin")
        call remap(dom, dom%topo%f_grnd,  gn, f_grnd_m,  gm, "bilin")
        call remap(dom, dom%topo%z_sl,    gn, z_sl_m,    gm, "bilin")
        call remap(dom, dom%topo%regions, gn, regions_m, gm, "nn")
        call remap(dom, dom%topo%basins,  gn, basins_m,  gm, "nn")

        call marshelf_update(dom%mshlf, H_ice_m, z_bed_m, f_grnd_m, regions_m, basins_m, &
                             z_sl_m, dx=dom%ctl%dx_mshlf)

        ! bmb_shlf / T_shlf now live on grid_mshlf in dom%mshlf%now;
        ! couple_marine_to_yelmo lands them on the Yelmo grid in the ice-sheet step.

    end subroutine step_marine_shelf_esm

    ! ---------------------------------------------------------------------------
    subroutine couple_esm_extras_to_yelmo(dom, esm)
        ! Land the ESM-owned Yelmo inputs that the shared couplers do not cover.
        ! Currently just subglacial discharge (Qd, Greenland frontal melt), remapped
        ! from the esm grid onto the Yelmo grid. Called before step_icesheet so Qd
        ! is in place when Yelmo runs, matching the paradigm used by every module.
        type(ice_domain),        intent(inout) :: dom
        type(esm_forcing_class), intent(inout) :: esm

        real(wp), allocatable :: Qd_y(:,:)
        character(len=256) :: ge, gy

        ge = trim(dom%ctl%grid_clim)    ! esm grid
        gy = trim(dom%ctl%grid_yelmo)

        call remap(dom, esm%Qd_ann, ge, Qd_y, gy, "con")
        dom%yelmo%bnd%Qd = Qd_y
    end subroutine couple_esm_extras_to_yelmo

end program yelmox_esm
