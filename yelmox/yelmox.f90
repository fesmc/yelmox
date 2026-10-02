program yelmox
    ! Multigrid yelmox driver (single domain).
    !
    ! Initializes one kryos_domain (each sub-model on its own configurable grid)
    ! plus the hi-res topography reference hub and the coupler maps, builds the
    ! initial boundary state (or restores a restart bundle), and runs the coupling
    ! time loop with per-module output. The multi-domain (bipolar) variant lives
    ! in yelmox_bipolar/. See docs/multigrid.md and libs/yelmox_domain.f90.

    use nml
    use timestepping
    use timeout
    use yelmo, only : yelmo_load_command_line_args, wp, yelmo_end
    use fastisostasy, only : bsl_class, bsl_init, bsl_update
    use kryos,          only : kryos_domain, domain_init
    use kryos_regions,  only : domain_regions_init
    use kryos_coupling, only : step_spinup_tuning, step_isostasy, couple_to_yelmo, &
                               step_icesheet, refresh_hub, step_climate, step_smb, &
                               step_marine_shelf
    use kryos_startup,  only : domain_startup, run_restart_write
    use kryos_forcing,  only : tsforcing_class, tsforcing_init, tsforcing_update, &
                               tsforcing_kill, tsforcing_restart_due, &
                               tsforcing_restart_fldr, tsforcing_restart_read, &
                               tsforcing_write_step
    use kryos_output,   only : domain_write_init, domain_write_step, &
                               domain_write_init_sm, domain_write_step_sm, &
                               domain_write_1D

    implicit none

    character(len=512) :: path_par
    type(tstep_class)  :: ts
    type(kryos_domain)   :: dom
    type(bsl_class)    :: bsl        ! shared, driver-owned barystatic sea level
    type(timeout_class) :: tm_2D, tm_2Dsm, tm_1D, tm_rst

    character(len=512) :: outfldr
    real(wp)           :: dtt

    ! Run control ([ctrl]): the group holding this run phase's timeline (e.g.
    ! "spinup", "transient"; "ctrl" = [ctrl] itself), and whether the timeline
    ! is in calendar years (tstep_const then a calendar time, against calendar_ref).
    character(len=56)  :: run_step
    logical            :: calendar
    real(wp)           :: calendar_ref

    ! Transient time-series forcing (tsgen), owned by the driver. The single
    ! forcing value f_now is mapped onto the snapclim anomalies via per-channel
    ! gains ([tsforcing]): dTa = f_now*f_ta, dTo = f_now*f_to, dSo = f_now*f_so.
    ! The tsforcing_class also owns the forcing-increment restart bookkeeping and
    ! the kill switch (see libs/yelmox_domain.f90).
    type(tsforcing_class) :: tsf
    real(wp) :: fvar

    ! Parameter file path from the command line (runme passes it per run).
    call yelmo_load_command_line_args(path_par)

    ! Timestepping (driver-owned; the [run_step] group holds the shared timeline).
    call nml_read(path_par, "ctrl", "run_step",     run_step)
    call nml_read(path_par, "ctrl", "calendar",     calendar)
    call nml_read(path_par, "ctrl", "calendar_ref", calendar_ref)
    call tstep_init(ts, path_par, trim(run_step), dtt, time_ref=calendar_ref, cal=calendar)

    ! Single-domain runs write to the run dir.
    outfldr = "./"

    ! Shared, driver-owned barystatic sea level (one per run).
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, ts%time_rel)

    ! Initialize the domain: sub-models + hi-res hub + coupler maps. The domain
    ! reads the timeline values it needs from the same [run_step] group.
    call domain_init(dom, path_par, ts%time, timeline_group=trim(run_step))

    ! Define regions of interest for 1D output (must precede the first yelmo_update).
    call domain_regions_init(dom, trim(outfldr))

    ! Transient time-series forcing (tsgen -> snapclim anomalies). Initialize
    ! before startup so the initial (cold-start) climate carries the same
    ! anomalies as the time loop. tsforcing reads [tsforcing] + [tsgen]; on a
    ! restart run, resume the series from the saved tsgen state in the bundle.
    call tsforcing_init(tsf, path_par, ts%time)
    if (trim(dom%ctl%restart) /= "None") call tsforcing_restart_read(tsf, trim(dom%ctl%restart))

    ! Cold start: build the initial boundary state. Restart: restore the bundle
    ! (incl. the shared bsl), rebuild the hi-res hub from the restored models,
    ! then re-establish the climate/smb and marine-shelf forcing from the
    ! restored state (the bundle does not hold them), so the first step and the
    ! first output see a valid boundary state.
    call domain_startup(dom, ts, bsl, tsf=tsf)
    if (trim(dom%ctl%restart) /= "None") then
        call step_climate(dom, ts, tsf)
        call step_smb(dom, ts)
        call step_marine_shelf(dom, ts)
    end if

    write(*,*)
    write(*,*) "yelmox: domain initialized"
    write(*,*) "  domain      : "//trim(dom%ctl%domain)
    write(*,*) "  Yelmo grid  : "//trim(dom%ctl%grid_ice), dom%yelmo%grd%G%nx, dom%yelmo%grd%G%ny
    write(*,*) "  topo grid   : "//trim(dom%ctl%grid_hub),  dom%topo%nx,      dom%topo%ny
    write(*,*) "  coupler maps: ", dom%cpl%nmaps
    write(*,*)

    ! === output + restart schedules (a restart bundle is always written at time_end) ===
    call timeout_init(tm_2D,   path_par, "tm_2D",   "heavy",   ts%time_init, ts%time_end)
    call timeout_init(tm_2Dsm, path_par, "tm_2Dsm", "medium",  ts%time_init, ts%time_end)
    call timeout_init(tm_1D,   path_par, "tm_1D",   "small",   ts%time_init, ts%time_end)
    call timeout_init(tm_rst,  path_par, "tm_rst",  "restart", ts%time_init, ts%time_end)

    ! Output files: 2D (one file per module, on its own grid), 2D small (reduced
    ! field set) and 1D timeseries.
    if (tm_2D%active)   call domain_write_init(dom, trim(outfldr), ts%time)
    if (tm_2Dsm%active) call domain_write_init_sm(dom, trim(outfldr), ts%time)
    if (tm_1D%active)   call domain_write_1D(dom, trim(outfldr), ts%time, init=.TRUE.)

    ! === main time loop ===
    ! Output and restarts are written at the top of the loop for the current
    ! time (time_init on the first pass), then the loop exits once the run is
    ! finished (time_end reached, or the kill switch tripped), else the time is
    ! advanced and the domain stepped. The final state and a final restart
    ! bundle are always written.
    call tstep_print_header(ts)
    do

        if (tm_2D%active .and. (timeout_check(tm_2D, ts%time) .or. ts%is_finished)) then
            call domain_write_step(dom, trim(outfldr), ts%time)
        end if

        if (tm_2Dsm%active .and. (timeout_check(tm_2Dsm, ts%time) .or. ts%is_finished)) then
            call domain_write_step_sm(dom, trim(outfldr), ts%time)
        end if

        if (tm_1D%active .and. (timeout_check(tm_1D, ts%time) .or. ts%is_finished)) then
            call domain_write_1D(dom, trim(outfldr), ts%time)
            call tsforcing_write_step(tsf, dom%yelmo%reg%fnm, ts%time)
        end if

        ! Restart bundle (domain + shared bsl + tsforcing state).
        if (timeout_check(tm_rst, ts%time) .or. ts%is_finished) then
            call run_restart_write(dom, bsl, ts%time, tsf=tsf)
        end if

        if (ts%is_finished) exit

        call tstep_update(ts, dtt)
        call tstep_print(ts)

        ! Shared sea level: update once per step, before the domain advances.
        call bsl_update(bsl, ts%time_rel)

        ! Transient forcing: advance the series every step (feedback methods need
        ! the response-derivative window); response variable = ice volume [Gt].
        if (tsf%active) then
            fvar = dom%yelmo%reg%V_ice * dom%yelmo%bnd%c%rho_ice * 1e-3_wp
            call tsforcing_update(tsf, ts%time, var=fvar)
        end if

        ! === coupling sequence ===
        call step_spinup_tuning(dom, ts)  ! relaxation ramp + cb_ref/tf_corr tuning (opt)
        call step_isostasy(dom, ts, bsl)  ! bedrock + sea level, this step
        call couple_to_yelmo(dom)         ! bedrock now; smb + shelf melt lag one step
        call step_icesheet(dom, ts)       ! yelmo_update
        call refresh_hub(dom)             ! hi-res geometry from the models
        call step_climate(dom, ts, tsf)   ! climate (dt_clim cadence)
        call step_smb(dom, ts)            ! surface mass balance
        call step_marine_shelf(dom, ts)   ! shelf melt

        ! Forcing-increment restart each |Δf| > restart_every_df (folders
        ! restart-<n>), so a ramp can be branched at fixed forcing levels.
        if (tsforcing_restart_due(tsf)) then
            write(*,*) "yelmox: forcing-increment restart at f =", tsf%tsg%f_now
            call run_restart_write(dom, bsl, ts%time, tsf=tsf, &
                                   fldr=trim(tsforcing_restart_fldr(tsf)))
        end if

        ! Kill switch: stop once the response has equilibrated at a forcing
        ! bound. The top of the loop then writes the final state and exits.
        if (tsforcing_kill(tsf)) then
            write(*,*) "yelmox: tsgen kill switch tripped at time =", ts%time
            ts%is_finished = .TRUE.
        end if
    end do

    write(*,*)
    write(*,*) "yelmox: run complete at time =", ts%time
    write(*,*) "  H_ice max   =", maxval(dom%yelmo%tpo%now%H_ice)

    ! Finalize Yelmo (deallocates model state) -- must come after the last
    ! access to dom%yelmo (it deallocates tpo%now%H_ice etc).
    call yelmo_end(dom%yelmo, time=ts%time)

end program yelmox
