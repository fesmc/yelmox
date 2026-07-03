program yelmox_mgesm
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
    ! The ESM forcing modules (esm + smbpal + marine_shelf) share one working
    ! grid, the "esm grid" = grid_mshlf: esm_clim_update / esm_forcing_update /
    ! esm_variability_update and the marine-shelf interpolation are co-gridded, as
    ! in the yelmox_esm.f90 monolith. Geometry is remapped from the hub onto that
    ! grid and the SMB / ocean boundary conditions are aggregated back to Yelmo.
    ! With grid_mshlf == grid_smb == grid_yelmo == grid_name every remap is an
    ! identity copy, reproducing yelmox_esm.f90.
    !
    ! Output (yelmo2D / yelmo2Dsm / yelmo1D_esm and the CMIP-formatted files) is
    ! kept identical to yelmox_esm.f90 via yelmox_esm_output. See docs/multigrid.md.

    use nml
    use ncio
    use timestepping
    use timeout
    use yelmo
    use ice_optimization

    use esm
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

    ! ESM run control (esm-specific config; the shared domain_ctl carries the rest)
    type esm_ctl_params
        character(len=56) :: run_step
        real(wp) :: time_init, time_end, dtt, time_equil
        real(wp) :: time_ref(2), time_hist(2), time_proj(2), time_esm_ref(2)
        character(len=56) :: clim_var
        integer  :: clim_seed
        character(len=56) :: tstep_method
        real(wp) :: tstep_const
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

    ! Run step selects which timeline / esm-timing group is read.
    call nml_read(path_par, "ctrl", "run_step", ec%run_step)

    ! [esm] group: experiment identity + physics parameters (into the esm object).
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

    ! [run_step] group: timeline + esm reference/history/projection periods.
    call nml_read(path_par, trim(ec%run_step), "time_init",    ec%time_init)
    call nml_read(path_par, trim(ec%run_step), "time_end",     ec%time_end)
    call nml_read(path_par, trim(ec%run_step), "dtt",          ec%dtt)
    call nml_read(path_par, trim(ec%run_step), "time_equil",   ec%time_equil)
    call nml_read(path_par, trim(ec%run_step), "time_ref",     ec%time_ref)
    call nml_read(path_par, trim(ec%run_step), "time_hist",    ec%time_hist)
    call nml_read(path_par, trim(ec%run_step), "time_proj",    ec%time_proj)
    call nml_read(path_par, trim(ec%run_step), "time_esm_ref", ec%time_esm_ref)
    call nml_read(path_par, trim(ec%run_step), "clim_var",     ec%clim_var)
    call nml_read(path_par, trim(ec%run_step), "clim_seed",    ec%clim_seed)
    call nml_read(path_par, trim(ec%run_step), "tstep_method", ec%tstep_method)
    call nml_read(path_par, trim(ec%run_step), "tstep_const",  ec%tstep_const)
    call nml_read(path_par, trim(ec%run_step), "kill_shelves", ec%kill_shelves)

    ! Seed the RNG for climate-variability reproducibility.
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = ec%clim_seed
    call random_seed(put=seed)

    ! Single-domain runs write to the run dir.
    outfldr     = "./"
    file2D      = trim(outfldr)//"yelmo2D.nc"
    file2Dsm    = trim(outfldr)//"yelmo2Dsm.nc"
    file1D_esm  = trim(outfldr)//"yelmo1D_esm.nc"
    file1D_cmip = trim(outfldr)//"yelmo1D_cmip.nc"
    file2D_cmip = trim(outfldr)//"yelmo2D_cmip.nc"
    file_isos   = trim(outfldr)//"fastisostasy.nc"
    file_bsl    = trim(outfldr)//"bsl.nc"

    write(*,*)
    write(*,*) "yelmox_mgesm: run_step = "//trim(ec%run_step)
    write(*,*) "  time_init/end/dtt: ", ec%time_init, ec%time_end, ec%dtt
    write(*,*) "  esm experiment:    "//trim(ec%experiment)//"  use_esm = ", ec%use_esm
    write(*,*)

    ! Timestepping (driver-owned).
    call tstep_init(ts, ec%time_init, ec%time_end, method=ec%tstep_method, units="year", &
                    time_ref=2000.0_wp, const_rel=0.0_wp, const_cal=ec%tstep_const)

    ! Shared, driver-owned barystatic sea level.
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, ts%time_rel)

    ! Initialize the domain (sub-models + hub + coupler maps), skipping snapclim:
    ! this driver supplies its own climate via the esm object.
    call domain_init(dom, path_par, ts%time, init_climate=.false.)

    ! Inject the driver-owned timeline values the domain logic needs.
    dom%ctl%tstep_method = ec%tstep_method
    dom%ctl%dtt          = ec%dtt
    dom%ctl%time_equil   = ec%time_equil

    ! Regions of interest for 1D output (must precede the first yelmo_update).
    call domain_regions_init(dom, trim(outfldr))

    ! Initialize the ESM forcing on the esm grid (= grid_mshlf, shared with mshlf).
    esm_path_par = trim(outfldr)//"/"//trim(ec%par_file)
    call esm_forcing_init(esm, esm_path_par, dom%ctl%domain, trim(dom%ctl%grid_mshlf), &
                          run_type=ec%run_step, gcm=ec%esm_name, experiment=ec%experiment, &
                          use_esm=ec%use_esm, use_smb=ec%use_smb, use_var=ec%use_var, &
                          use_hist=ec%use_hist, use_proj=ec%use_proj)

    write(*,*)
    write(*,*) "yelmox_mgesm: domain initialized"
    write(*,*) "  domain      : "//trim(dom%ctl%domain)
    write(*,*) "  Yelmo grid  : "//trim(dom%ctl%grid_yelmo), dom%yelmo%grd%nx, dom%yelmo%grd%ny
    write(*,*) "  esm  grid   : "//trim(dom%ctl%grid_mshlf)
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
        call bsl_restart_read(bsl, trim(dom%ctl%restart)//"/bsl_restart.nc")
        call bsl_update(bsl, ts%time_rel)
        call domain_restart_read(dom, trim(dom%ctl%restart), ts, bsl)
        call refresh_htopo(dom)
        call step_climate_esm(dom, esm, ec, ts)
        call step_marine_shelf_esm(dom, esm, ec, ts)
    end if

    ! ================= OUTPUT SETUP ==========================================

    call timeout_init(tm_1D,   path_par, "tm_1D",   "small",  ec%time_init, ec%time_end)
    call timeout_init(tm_2D,   path_par, "tm_2D",   "heavy",  ec%time_init, ec%time_end)
    call timeout_init(tm_2Dsm, path_par, "tm_2Dsm", "medium", ec%time_init, ec%time_end)

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
        call step_optimize(dom, ts)
        call step_isostasy(dom, ts, bsl)
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
            if (mod(nint(ts%time_elapsed*100), nint(ec%dt_formatted*100)) == 0) then
                call write_step_2D_cmip(dom%yelmo, dom%mshlf, file2D_cmip, ts%time)
                call write_step_1D_cmip(dom%yelmo, dom%mshlf, file1D_cmip, ts%time)
            end if
        end if

        ! === Restart bundle ===
        if (dom%ctl%dt_restart > 0.0_wp .and. &
            mod(nint(ts%time*100), nint(dom%ctl%dt_restart*100)) == 0) then
            call domain_restart_write(dom, ts%time)
            call restart_bundle_mkdir(ts%time)
            call bsl_restart_write(bsl, trim(restart_bundle_dir(ts%time))//"/bsl_restart.nc", ts%time)
        end if

    end do

    ! Final state + restart bundle (incl. shared bsl).
    call write_step_2D_combined(dom%yelmo, dom%isos, esm, dom%mshlf, dom%smb, &
                                ec%use_smb, file2D, ts%time)
    call yelmo_regions_write(dom%yelmo, ts%time)
    call domain_restart_write(dom, ts%time)
    call restart_bundle_mkdir(ts%time)
    call bsl_restart_write(bsl, trim(restart_bundle_dir(ts%time))//"/bsl_restart.nc", ts%time)

    write(*,*)
    write(*,*) "yelmox_mgesm: run complete at time =", ts%time
    write(*,*) "  H_ice max =", maxval(dom%yelmo%tpo%now%H_ice)

    call yelmo_end(dom%yelmo, time=ts%time)
    deallocate(seed)

contains

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
        real(wp), allocatable :: z_bed_i(:,:), H_ice_i(:,:), z_bed_y(:,:), z_ss_y(:,:)
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
        call remap(dom, dom%isos%out%z_bed, gi, z_bed_y, gy, "con")
        call remap(dom, dom%isos%out%z_ss,  gi, z_ss_y,  gy, "con")
        dom%yelmo%bnd%z_bed = z_bed_y
        dom%yelmo%bnd%z_sl  = z_ss_y

        ! Hub mirror + first ESM forcing (init=.true. runs the smbpal ITM equil).
        call refresh_htopo(dom)
        call step_climate_esm(dom, esm, ec, ts, init=.true.)
        call step_marine_shelf_esm(dom, esm, ec, ts)

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
        ! ESM climate + surface mass balance. The esm forcing runs on the esm grid
        ! (grid_mshlf): geometry is remapped from the hub, the three esm updates
        ! produce the atmospheric/ocean anomaly fields, smbpal (or the direct-SMB
        ! path) is evaluated on grid_smb, and smb / T_srf are aggregated back to the
        ! Yelmo grid. The ocean anomalies (esm%dto/dso, esm%to_ref/so_ref) are left
        ! for step_marine_shelf_esm, which shares the esm grid.
        type(ice_domain),        intent(inout) :: dom
        type(esm_forcing_class), intent(inout) :: esm
        type(esm_ctl_params),    intent(in)    :: ec
        type(tstep_class),       intent(in)    :: ts
        logical, intent(in), optional          :: init

        real(wp), allocatable :: z_srf_f(:,:), H_ice_f(:,:), z_bed_f(:,:)
        real(wp), allocatable :: f_grnd_f(:,:), z_sl_f(:,:), basins_f(:,:)
        real(wp), allocatable :: tas_s(:,:,:), pr_s(:,:,:), z_srf_s(:,:), H_ice_s(:,:)
        real(wp), allocatable :: smb_y(:,:), tsrf_y(:,:), Qd_y(:,:)
        character(len=256) :: gf, gn, gy, gs
        logical :: is_init

        is_init = .false.
        if (present(init)) is_init = init

        gf = trim(dom%ctl%grid_mshlf)   ! esm grid (shared with marine_shelf)
        gn = trim(dom%ctl%grid_name)    ! hub grid
        gy = trim(dom%ctl%grid_yelmo)
        gs = trim(dom%ctl%grid_smb)

        ! Geometry: hub -> esm grid.
        call remap(dom, dom%topo%z_srf,   gn, z_srf_f,  gf, "bilin")
        call remap(dom, dom%topo%H_ice,   gn, H_ice_f,  gf, "bilin")
        call remap(dom, dom%topo%z_bed,   gn, z_bed_f,  gf, "bilin")
        call remap(dom, dom%topo%f_grnd,  gn, f_grnd_f, gf, "bilin")
        call remap(dom, dom%topo%z_sl,    gn, z_sl_f,   gf, "bilin")
        call remap(dom, dom%topo%basins,  gn, basins_f, gf, "nn")

        ! Step 1: reference climatology (lapse-rate / precip scaling to z_srf).
        call esm_clim_update(esm, z_srf_f, ts%time, ec%time_ref, ec%use_smb, &
                             dom%ctl%domain, gf)

        ! Extrapolate reference ocean into ice-shelf interiors.
        if (dom%mshlf%par%extrap_shlf) then
            call ocn_variable_extrapolation(esm%to_ref%var(:,:,:,1), H_ice_f, basins_f, &
                                            -esm%to_var_ref%z, z_bed_f)
            call ocn_variable_extrapolation(esm%so_ref%var(:,:,:,1), H_ice_f, basins_f, &
                                            -esm%so_var_ref%z, z_bed_f)
        end if

        ! Step 2: anomaly fields (historical / projection / homogeneous).
        call esm_forcing_update(esm, dom%mshlf, ts%time, ec%use_esm, ec%time_ref, &
                                ec%time_hist, ec%time_proj, ec%time_esm_ref, dom%ctl%domain, &
                                H_ice_f, basins_f, z_bed_f, f_grnd_f, z_sl_f, ec%use_smb, &
                                use_ref_atm=.false., use_ref_ocn=.false.)

        ! Step 3: variability anomaly.
        call esm_variability_update(esm, dom%mshlf, ts%time, ec%dtt, ec%clim_var, ec%time_ref, &
                                    H_ice_f, basins_f, z_bed_f, f_grnd_f, z_sl_f, ec%use_var, &
                                    use_ref_atm=.false., use_ref_ocn=.false.)

        ! === Atmospheric boundary conditions (SMB on grid_smb) ===
        if (ec%use_smb) then
            ! Direct SMB from the esm fields (co-gridded on the esm grid gf==gs
            ! for the standard ESM setup); assign into the smbpal annual fields.
            dom%smb%ann%smb  = esm%smb_ann + sum(esm%dsmb, dim=3)/12.0_wp &
                               - esm%dsmbdz*(dom%yelmo%dta%pd%z_srf - dom%yelmo%tpo%now%z_srf)
            dom%smb%ann%tsrf = sum(esm%t2m + esm%dts + esm%dts_var, dim=3)/12.0_wp
            where(dom%yelmo%tpo%now%H_ice > 0.0_wp .and. &
                  sum(esm%t2m + esm%dts + esm%dts_var, dim=3)/12.0_wp > 273.15_wp) &
                dom%smb%ann%tsrf = 273.15_wp
        else
            ! smbpal on grid_smb: atmospheric forcing from the esm grid, geometry
            ! from the hub. init=.true. runs the ITM snowpack equilibration first.
            call remap(dom, esm%t2m + esm%dts + esm%dts_var, gf, tas_s, gs, "bilin")
            call remap(dom, esm%pr  * esm%dpr * esm%dpr_var, gf, pr_s,  gs, "bilin")
            call remap(dom, dom%topo%z_srf, gn, z_srf_s, gs, "bilin")
            call remap(dom, dom%topo%H_ice, gn, H_ice_s, gs, "bilin")
            if (is_init .and. trim(dom%smb%par%abl_method) == "itm") then
                call remap(dom, esm%t2m + esm%dts, gf, tas_s, gs, "bilin")
                call remap(dom, esm%pr  * esm%dpr, gf, pr_s,  gs, "bilin")
                call smbpal_update_monthly_equil(dom%smb, tas_s, pr_s, z_srf_s, H_ice_s, &
                                                 ts%time_rel, time_equil=100.0_wp)
                call remap(dom, esm%t2m + esm%dts + esm%dts_var, gf, tas_s, gs, "bilin")
                call remap(dom, esm%pr  * esm%dpr * esm%dpr_var, gf, pr_s,  gs, "bilin")
            end if
            call smbpal_update_monthly(dom%smb, tas_s, pr_s, z_srf_s, H_ice_s, ts%time_rel)
        end if

        ! Aggregate SMB / surface temperature -> Yelmo grid (conservative).
        call remap(dom, dom%smb%ann%smb,  gs, smb_y,  gy, "con")
        call remap(dom, dom%smb%ann%tsrf, gs, tsrf_y, gy, "con")
        dom%yelmo%bnd%smb   = smb_y * dom%yelmo%bnd%c%conv_we_ie * 1e-3
        dom%yelmo%bnd%T_srf = tsrf_y

        ! Subglacial discharge (Greenland frontal melt) -> Yelmo grid.
        call remap(dom, esm%Qd_ann, gf, Qd_y, gy, "con")
        dom%yelmo%bnd%Qd = Qd_y

    end subroutine step_climate_esm

    ! ---------------------------------------------------------------------------
    subroutine step_marine_shelf_esm(dom, esm, ec, ts)
        ! ESM ocean boundary conditions + marine-shelf basal melt on the esm grid
        ! (grid_mshlf): interpolate the reference ocean to shelf depth, add the esm
        ! ocean anomalies (set by step_climate_esm), run marshelf_update, and
        ! aggregate bmb_shlf / T_shlf back to the Yelmo grid.
        type(ice_domain),        intent(inout) :: dom
        type(esm_forcing_class), intent(inout) :: esm
        type(esm_ctl_params),    intent(in)    :: ec
        type(tstep_class),       intent(in)    :: ts

        real(wp), allocatable :: H_ice_f(:,:), z_bed_f(:,:), f_grnd_f(:,:), z_sl_f(:,:)
        real(wp), allocatable :: regions_f(:,:), basins_f(:,:)
        real(wp), allocatable :: bmb_y(:,:), Tshlf_y(:,:)
        character(len=256) :: gf, gn, gy

        if (.not. dom%ctl%with_marine_shelf) return

        gf = trim(dom%ctl%grid_mshlf)
        gn = trim(dom%ctl%grid_name)
        gy = trim(dom%ctl%grid_yelmo)

        ! Geometry + masks: hub -> esm grid.
        call remap(dom, dom%topo%H_ice,   gn, H_ice_f,   gf, "bilin")
        call remap(dom, dom%topo%z_bed,   gn, z_bed_f,   gf, "bilin")
        call remap(dom, dom%topo%f_grnd,  gn, f_grnd_f,  gf, "bilin")
        call remap(dom, dom%topo%z_sl,    gn, z_sl_f,    gf, "bilin")
        call remap(dom, dom%topo%regions, gn, regions_f, gf, "nn")
        call remap(dom, dom%topo%basins,  gn, basins_f,  gf, "nn")

        ! Reference ocean interpolated to shelf depth, plus esm anomalies.
        call marshelf_interp_shelf(dom%mshlf%now%T_shlf, dom%mshlf, esm%to_ref%var(:,:,:,1), &
                                   H_ice_f, z_bed_f, f_grnd_f, z_sl_f, -esm%to_ref%z)
        call marshelf_interp_shelf(dom%mshlf%now%S_shlf, dom%mshlf, esm%so_ref%var(:,:,:,1), &
                                   H_ice_f, z_bed_f, f_grnd_f, z_sl_f, -esm%so_ref%z)
        dom%mshlf%now%T_shlf = dom%mshlf%now%T_shlf + esm%dto + esm%dto_var
        dom%mshlf%now%S_shlf = dom%mshlf%now%S_shlf + esm%dso + esm%dso_var

        if (trim(dom%ctl%domain) == "Greenland") then
            dom%mshlf%now%dT_shlf   = dom%mshlf%now%T_shlf + esm%dto
            dom%mshlf%par%tf_method = 2
        end if

        call marshelf_update(dom%mshlf, H_ice_f, z_bed_f, f_grnd_f, regions_f, basins_f, &
                             z_sl_f, dx=dom%ctl%dx_mshlf)

        ! Aggregate outputs -> Yelmo grid (conservative).
        call remap(dom, dom%mshlf%now%bmb_shlf, gf, bmb_y,   gy, "con")
        call remap(dom, dom%mshlf%now%T_shlf,   gf, Tshlf_y, gy, "con")
        dom%yelmo%bnd%bmb_shlf = bmb_y
        dom%yelmo%bnd%T_shlf   = Tshlf_y

    end subroutine step_marine_shelf_esm

end program yelmox_mgesm
