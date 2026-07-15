program yelmox_rembo
    ! REMBO-coupled yelmox driver (Greenland, single grid).
    !
    ! Reuses the multigrid ice_domain (yelmox_domain) for the ice sheet, isostasy,
    ! ocean (snapclim), marine_shelf, sediments, geothermal, the hi-res hub and the
    ! shared barystatic sea level, and drives them with the shared step_* / coupler
    ! primitives. REMBO replaces the snapclim+smbpal atmosphere/SMB: step_rembo runs
    ! REMBO on the domain grid and stages its smb/T_srf into the SMB carrier
    ! (dom%smb%ann), which couple_smb_to_yelmo lands on the Yelmo grid like any other
    ! module. Hysteresis forcing (dT) and REMBO itself are driver-owned (REMBO keeps
    ! its module-global state rembo_ann); see yelmox_rembo/legacy for the original.

    use nml
    use ncio
    use timestepping
    use timeout
    use yelmo, only : yelmo_load_command_line_args, wp, dp, yelmo_end, &
                      yelmo_init_state, yelmo_update_equil, yelmo_print_bound
    use fastisostasy, only : bsl_class, bsl_init, bsl_update, &
                             isos_init_ref, isos_init_state
    use yelmox_climate, only : climate_update
    use rembo_sclimate, only : rembo_init, rembo_update, rembo_equilibrate, &
                               rembo_ann, rembo_restart_write
    use tsgen, only : tsgen_class
    use yelmox_domain
    use yelmox_rembo_output

    implicit none

    character(len=512)  :: path_par
    type(tstep_class)   :: ts
    type(ice_domain)    :: dom
    type(bsl_class)     :: bsl          ! shared, driver-owned barystatic sea level
    type(tsforcing_class) :: tsf        ! driver-owned transient forcing (tsgen)
    type(timeout_class) :: tm_2D, tm_2Dsm, tm_1D

    character(len=512)  :: outfldr, file_rembo2D, file_rembo1D
    real(wp) :: time_equil, dtt, dtt_now, deltat_tot
    logical  :: write_restart
    real(wp) :: dT_summer, dT_ann, dT_ocn

    ! Parameter file path from the command line (runme passes it per run).
    call yelmo_load_command_line_args(path_par)

    ! --- run control ([ctrl]: shared timeline + REMBO-flavor switches) ---
    call timeline_init(ts, dtt, path_par, "ctrl")
    call nml_read(path_par, "ctrl", "time_equil",   time_equil)
    call nml_read(path_par, "ctrl", "write_restart", write_restart)

    outfldr      = "./"
    file_rembo2D = trim(outfldr)//"rembo.nc"
    file_rembo1D = trim(outfldr)//"rembo_ts.nc"

    ! Shared, driver-owned barystatic sea level (one per run).
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, ts%time_rel)

    ! Domain: ice sheet + isostasy + ocean (snapclim) + smbpal carrier + marine
    ! shelf + sediments + geothermal + hi-res hub + coupler maps. The domain
    ! reads the timeline values it needs from the same [ctrl] group.
    call domain_init(dom, path_par, ts%time)

    ! REMBO climate + hysteresis (driver-owned; not part of ice_domain). REMBO keeps
    ! its module-global state in rembo_ann and loads its own parameters (options_rembo).
    call rembo_init(real(ts%time, dp))
    call tsforcing_init(tsf, path_par, ts%time)
    if (trim(dom%ctl%restart) /= "None") call tsforcing_restart_read(tsf, trim(dom%ctl%restart))

    ! Regions of interest for 1D output (must precede the first yelmo_update).
    call domain_regions_init(dom, trim(outfldr))

    ! Initial transient forcing (dT_summer / dT_ann / dT_ocn).
    call update_forcing(tsf, dom, ts, dT_summer, dT_ann, dT_ocn)

    ! Cold start: build the initial boundary state. Restart: restore the bundle
    ! (incl. the shared bsl) and re-establish the REMBO + ocean forcing.
    if (trim(dom%ctl%restart) == "None") then
        call rembo_cold_start(dom, ts, bsl, tsf, dT_summer, dT_ocn, time_equil, dtt)
    else
        call domain_startup(dom, ts, bsl)
        call step_rembo(dom, ts, tsf, dT_summer, dT_ocn)
        call step_marine_shelf(dom, ts)
    end if

    write(*,*)
    write(*,*) "yelmox_rembo: domain initialized"
    write(*,*) "  domain      : "//trim(dom%ctl%domain)
    write(*,*) "  Yelmo grid  : "//trim(dom%ctl%grid_yelmo), dom%yelmo%grd%G%nx, dom%yelmo%grd%G%ny
    write(*,*)

    ! === output setup ===
    ! Standard per-module files (yelmo.nc / yelmo_sm.nc / yelmo_ts.nc / isos.nc /
    ! mshlf.nc / htopo.nc) via the shared routines, exactly as the other flavors;
    ! REMBO-specific fields go to rembo.nc (2D) and rembo_ts.nc (1D).
    call timeout_init(tm_2D,   path_par, "tm_2D",   "heavy",  ts%time_init, ts%time_end)
    call timeout_init(tm_2Dsm, path_par, "tm_2Dsm", "medium", ts%time_init, ts%time_end)
    call timeout_init(tm_1D,   path_par, "tm_1D",   "small",  ts%time_init, ts%time_end)

    if (tm_2D%active) then
        call domain_write_init(dom, trim(outfldr), ts%time)
        call domain_write_step(dom, trim(outfldr), ts%time)
        call rembo_write_2D_init(dom%yelmo, file_rembo2D, ts%time, "years")
        call rembo_write_2D_step(rembo_ann, file_rembo2D, ts%time)
    end if
    if (tm_2Dsm%active) then
        call domain_write_init_sm(dom, trim(outfldr), ts%time)
        call domain_write_step_sm(dom, trim(outfldr), ts%time)
    end if
    if (tm_1D%active) then
        call domain_write_1D(dom, trim(outfldr), ts%time, init=.TRUE.)
        call rembo_write_1D_init(file_rembo1D, ts%time, "years", tsf%tsg)
        call rembo_write_1D_step(dom%yelmo, tsf%tsg, rembo_ann, file_rembo1D, ts%time, dT_ann, dT_ocn)
    end if

    ! Initial restart bundle.
    if (write_restart) call write_rembo_restart(dom, bsl, ts, tsf)

    ! === main time loop ===
    call tstep_print_header(ts)
    do while (.not. ts%is_finished)

        ! Transient experiments: shrink dtt (and REMBO's emb cadence) during the
        ! hysteresis ramp, restore afterwards.
        dtt_now = dtt
        if (tsf%active) then
            select case(trim(tsf%tsg%par%method))
                case("ramp-time","ramp-time-step")
                    deltat_tot = tsf%tsg%par%dt_init + tsf%tsg%par%dt_ramp + tsf%tsg%par%dt_conv + 100.0_wp
                    if (ts%time_elapsed < deltat_tot) then
                        dtt_now = min(5.0_wp, dtt)
                        rembo_ann%par%dtime_emb = real(dtt_now, dp)
                    else
                        dtt_now = dtt
                        rembo_ann%par%dtime_emb = 100.0_dp
                    end if
            end select
        end if

        call tstep_update(ts, dtt_now)
        call tstep_print(ts)

        ! Transient forcing for this step (from the tsgen series).
        call update_forcing(tsf, dom, ts, dT_summer, dT_ann, dT_ocn)

        ! Shared sea level: update once per step, before the domain advances.
        call bsl_update(bsl, ts%time_rel)

        ! Coupling sequence (REMBO replaces the generic climate step).
        call step_optimize(dom, ts)
        call step_isostasy(dom, ts, bsl)
        call step_icesheet(dom, ts)      ! couplers (smb/isos/marine) + yelmo_update
        call refresh_htopo(dom)
        call step_rembo(dom, ts, tsf, dT_summer, dT_ocn)  ! REMBO atmosphere/smb + snapclim ocean
        call step_marine_shelf(dom, ts)

        ! === output ===
        if (tm_2D%active .and. timeout_check(tm_2D, ts%time)) then
            call domain_write_step(dom, trim(outfldr), ts%time)
            call rembo_write_2D_step(rembo_ann, file_rembo2D, ts%time)
        end if
        if (tm_2Dsm%active .and. timeout_check(tm_2Dsm, ts%time)) &
            call domain_write_step_sm(dom, trim(outfldr), ts%time)
        if (tm_1D%active .and. timeout_check(tm_1D, ts%time)) then
            call domain_write_1D(dom, trim(outfldr), ts%time)
            call rembo_write_1D_step(dom%yelmo, tsf%tsg, rembo_ann, file_rembo1D, ts%time, dT_ann, dT_ocn)
        end if

        ! === restart bundle ===
        if (write_restart .and. tstep_due(ts%time, dom%ctl%dt_restart)) then
            call write_rembo_restart(dom, bsl, ts, tsf)
        end if

        ! tsgen kill switch (response equilibrated at a forcing bound).
        if (tsforcing_kill(tsf)) then
            write(*,"(a,f12.3,a,f12.3)") "tsgen:: kill switch activated. [time, f_now] = ", &
                    ts%time, ", ", tsf%tsg%f_now
            exit
        end if
    end do

    ! Final state + restart bundle.
    if (tm_2D%active) then
        call domain_write_step(dom, trim(outfldr), ts%time)
        call rembo_write_2D_step(rembo_ann, file_rembo2D, ts%time)
    end if
    if (tm_2Dsm%active) call domain_write_step_sm(dom, trim(outfldr), ts%time)
    if (tm_1D%active) then
        call domain_write_1D(dom, trim(outfldr), ts%time)
        call rembo_write_1D_step(dom%yelmo, tsf%tsg, rembo_ann, file_rembo1D, ts%time, dT_ann, dT_ocn)
    end if
    if (write_restart) call write_rembo_restart(dom, bsl, ts, tsf)

    write(*,*)
    write(*,*) "yelmox_rembo: run complete at time =", ts%time
    write(*,*) "  H_ice max   =", maxval(dom%yelmo%tpo%now%H_ice)

    call yelmo_end(dom%yelmo, time=ts%time)

contains

    subroutine update_forcing(tsf, dom, ts, dT_summer, dT_ann, dT_ocn)
        ! Update the transient temperature anomalies from the tsgen forcing series.
        ! REMBO uses its own channel mapping (not the snapclim dTa/dTo/dSo): f_now
        ! drives dT_summer, and dT_ann uses the REMBO winter factor (T_wintfac=1.6):
        !   dT_ann = 0.5*((1.6)*dT_summer + (1.0)*dT_summer) = 1.3*dT_summer.
        type(tsforcing_class), intent(inout) :: tsf
        type(ice_domain),      intent(in)    :: dom
        type(tstep_class),     intent(in)    :: ts
        real(wp),              intent(out)   :: dT_summer, dT_ann, dT_ocn

        real(wp) :: var

        if (tsf%active) then
            var = dom%yelmo%reg%V_ice * dom%yelmo%bnd%c%rho_ice * 1e-3_wp
            call tsforcing_update(tsf, ts%time, var=var)
            dT_summer = tsf%tsg%f_now * tsf%f_ta
            dT_ann    = 1.3_wp * dT_summer
            dT_ocn    = dT_ann * tsf%f_to
        else
            dT_summer = 0.0_wp
            dT_ann    = 0.0_wp
            dT_ocn    = 0.0_wp
        end if
    end subroutine update_forcing

    subroutine step_rembo(dom, ts, tsf, dT_summer, dT_ocn, init)
        ! Advance REMBO one coupling step: geometry from the hub, REMBO atmosphere +
        ! surface mass balance staged into the SMB carrier (dom%smb%ann), and ocean
        ! forcing from snapclim (with the optional hysteresis ocean anomaly). REMBO
        ! runs on grid_clim (== the Yelmo grid for the single-grid Greenland setup);
        ! couple_smb_to_yelmo applies the we->ie scaling + lim_pd_ice, on the Yelmo
        ! grid, inside step_icesheet. init=.true. equilibrates REMBO before the first
        ! update (cold start only).
        type(ice_domain),      intent(inout) :: dom
        type(tstep_class),     intent(in)    :: ts
        type(tsforcing_class), intent(in)    :: tsf
        real(wp),              intent(in)    :: dT_summer, dT_ocn
        logical, intent(in), optional :: init

        real(wp), allocatable :: z_srf_c(:,:), H_ice_c(:,:), z_sl_c(:,:), basins_c(:,:)
        character(len=256) :: gc, gn
        logical :: is_init

        is_init = .false.
        if (present(init)) is_init = init

        gc = trim(dom%ctl%grid_clim)
        gn = trim(dom%ctl%grid_name)

        ! Geometry from the hub -> REMBO/clim grid.
        call remap(dom, dom%topo%z_srf,  gn, z_srf_c,  gc, "bilin")
        call remap(dom, dom%topo%H_ice,  gn, H_ice_c,  gc, "bilin")
        call remap(dom, dom%topo%z_sl,   gn, z_sl_c,   gc, "bilin")
        call remap(dom, dom%topo%basins, gn, basins_c, gc, "nn")

        ! REMBO atmosphere + surface mass balance (double precision internally).
        if (is_init .and. .not. dom%yelmo%par%use_restart) then
            call rembo_equilibrate(real(ts%time, dp), real(z_srf_c, dp), real(H_ice_c, dp), &
                                   real(z_sl_c, dp), time_tot=10.0_dp)
        else
            call rembo_update(real(ts%time, dp), real(ts%time_rel, dp), real(dT_summer, dp), &
                              real(z_srf_c, dp), real(H_ice_c, dp), real(z_sl_c, dp))
        end if

        ! Stage the REMBO output into the SMB carrier (grid_smb == grid_clim here).
        dom%smb%ann%smb  = real(rembo_ann%smb,   wp)
        dom%smb%ann%tsrf = real(rembo_ann%T_srf, wp)

        ! Ocean forcing via the climate backend (grid_clim); optional hysteresis anomaly.
        call climate_update(dom%cl, dom%clim, z_srf=z_srf_c, time=ts%time, &
                            domain=dom%ctl%domain, dx=dom%ctl%dx_clim, basins=basins_c)
        if (tsf%active .and. trim(dom%cl%snp%par%ocn_type) == "const") &
            dom%clim%now%to_ann = dom%clim%now%to_ann + dT_ocn
    end subroutine step_rembo

    subroutine rembo_cold_start(dom, ts, bsl, tsf, dT_summer, dT_ocn, time_equil, dtt)
        ! Build the initial boundary state for a cold start (no restart bundle):
        ! isostasy reference/state, the first REMBO + ocean + marine forcing, the
        ! Yelmo state init, and the thermodynamic/dynamic equilibration passes.
        type(ice_domain),      intent(inout) :: dom
        type(tstep_class),     intent(in)    :: ts
        type(bsl_class),       intent(inout) :: bsl
        type(tsforcing_class), intent(in)    :: tsf
        real(wp),              intent(in)    :: dT_summer, dT_ocn
        real(wp),              intent(in)    :: time_equil, dtt

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

        ! Hub mirror + first REMBO / ocean / marine forcing.
        call refresh_htopo(dom)
        call step_rembo(dom, ts, tsf, dT_summer, dT_ocn, init=.true.)
        call step_marine_shelf(dom, ts)

        ! Assemble the Yelmo boundary state from the freshly produced outputs.
        call couple_smb_to_yelmo(dom)
        call couple_marine_to_yelmo(dom)

        ! Basal-friction optimization cold-start guess (equil_method == "opt").
        if (trim(dom%ctl%equil_method) == "opt") then
            dom%yelmo%dyn%par%till_method = -1
            if (.not. dom%yelmo%par%use_restart) then
                if (dom%opt%cf_init > 0.0_wp) then
                    dom%yelmo%dyn%now%cb_ref = dom%opt%cf_init
                else
                    dom%yelmo%dyn%now%cb_ref = dom%yelmo%dyn%now%cb_tgt
                end if
            end if
        end if

        ! Initialize Yelmo state variables (cold base).
        call yelmo_print_bound(dom%yelmo%bnd)
        call yelmo_init_state(dom%yelmo, time=ts%time, thrm_method="robin-cold")

        ! Optional LGM-like marine ice at the start.
        if (dom%ctl%greenland_init_marine_H) &
            dom%yelmo%tpo%now%H_ice = dom%yelmo%tpo%now%H_ice * 1.2_wp

        ! Equilibrate thermodynamics/dynamics (cold start only).
        if (.not. dom%yelmo%par%use_restart .and. dom%ctl%with_ice_sheet) then
            call yelmo_update_equil(dom%yelmo, ts%time, time_tot=10.0_wp, dt=1.0_wp, topo_fixed=.FALSE.)
            call yelmo_update_equil(dom%yelmo, ts%time, time_tot=time_equil, dt=dtt, topo_fixed=.TRUE.)
        end if
    end subroutine rembo_cold_start

    subroutine write_rembo_restart(dom, bsl, ts, tsf)
        ! Restart bundle: the shared domain sub-models + shared bsl
        ! (run_restart_write) plus REMBO's own restart, all in one folder.
        type(ice_domain),      intent(inout) :: dom
        type(bsl_class),       intent(inout) :: bsl
        type(tstep_class),     intent(in)    :: ts
        type(tsforcing_class), intent(inout) :: tsf

        call run_restart_write(dom, bsl, ts%time, tsf=tsf)
        call rembo_restart_write(trim(restart_bundle_dir(ts%time))//"/rembo_restart.nc", &
                real(ts%time, dp), real(dom%yelmo%tpo%now%z_srf, dp), &
                real(dom%yelmo%tpo%now%H_ice, dp), real(dom%yelmo%bnd%z_sl, dp))
    end subroutine write_rembo_restart

end program yelmox_rembo
