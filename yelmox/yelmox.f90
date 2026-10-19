program yelmox
    ! Multigrid yelmox driver (single domain).
    !
    ! Initializes one kryos_domain (each sub-model on its own configurable grid)
    ! plus the hi-res topography reference hub and the coupler maps, builds the
    ! initial boundary state (or restores a restart bundle), and runs the coupling
    ! time loop with per-module output. The multi-domain (bipolar) variant lives
    ! in yelmox_bipolar/. See docs/multigrid.md and libs/kryos*.f90.

    use nml
    use timestepping
    use timeout
    use precision, only : dp
    use yelmo, only : yelmo_load_command_line_args, wp, yelmo_end
    use fastisostasy, only : bsl_class, bsl_init, bsl_update
    use kryos,          only : kryos_domain, domain_init
    use yelmox_climate, only : climate_end
    use kryos_regions,  only : domain_regions_init
    use kryos_coupling, only : step_relax, step_optimize, step_isostasy, couple_to_yelmo, &
                               step_ice_bnd, step_icesheet, couple_yelmo_to_htopo, step_climate, &
                               step_surface, step_shelf
    use kryos_startup,  only : domain_startup, run_restart_write
    use kryos_forcing,  only : indices_class, indices_init, indices_update, indices_kill, &
                               indices_restart_due, indices_restart_fldr, &
                               indices_restart_read, indices_write_step
    use kryos_output,   only : domain_write_init, domain_write_step, &
                               domain_write_init_sm, domain_write_step_sm, &
                               domain_write_1D, domain_write_cmip
    use timer,          only : timer_class, timer_init, timer_step, &
                               timer_print_totals, timer_write_totals

    implicit none

    character(len=512) :: path_par
    type(tstep_class)  :: ts
    type(kryos_domain)   :: dom
    type(bsl_class)    :: bsl        ! shared, driver-owned barystatic sea level
    type(timeout_class) :: tm_2D, tm_2Dsm, tm_1D, tm_rst

    character(len=512) :: outfldr
    real(dp)           :: dtt

    ! Run control ([ctrl]): the group holding this run phase's timeline (e.g.
    ! "spinup", "transient"; "ctrl" = [ctrl] itself), and whether the timeline
    ! is in calendar years (tstep_const then a calendar time, against calendar_ref).
    character(len=56)  :: run_step
    logical            :: calendar
    real(dp)           :: calendar_ref

    ! Restart bundle folder to start from ([ctrl]); "None" = cold start.
    character(len=512) :: restart

    ! Indices ([indices], one tsgen series [idx_<name>] each), owned by the
    ! driver and advanced once per step on the forcing clock time_rel; read by
    ! name by the climate (snapshots axes, the modifiers such as [anom]). The
    ! registry also owns the kill switch and the forcing-increment restarts
    ! (see libs/kryos_forcing.f90).
    type(indices_class)   :: idx
    real(wp) :: fvar

    ! Wall time of each part of the time loop (fesm-utils timer), accumulated
    ! over the run: printed at the end, and written to timing.txt (tm_1D cadence).
    type(timer_class) :: tmr
    integer, parameter :: tt_output = 1, tt_forcing = 2, tt_relax = 3, tt_optimize = 4, &
                          tt_isostasy = 5, tt_couple_yelmo = 6, tt_icesheet = 7, &
                          tt_couple_htopo = 8, tt_climate = 9, tt_surface = 10, tt_shelf = 11
    logical :: timing_init = .TRUE.   ! timing.txt is created on the first 1D write

    ! Parameter file path from the command line (runme passes it per run).
    call yelmo_load_command_line_args(path_par)

    ! Timestepping (driver-owned; the [run_step] group holds the shared timeline).
    call nml_read(path_par, "ctrl", "run_step",     run_step)
    call nml_read(path_par, "ctrl", "calendar",     calendar)
    call nml_read(path_par, "ctrl", "calendar_ref", calendar_ref)
    call nml_read(path_par, "ctrl", "restart",      restart)
    call tstep_init(ts, path_par, trim(run_step), dtt, time_ref=calendar_ref, cal=calendar)

    ! Single-domain runs write to the run dir.
    outfldr = "./"

    ! Shared, driver-owned barystatic sea level (one per run).
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, ts%time_rel)

    ! Indices on the forcing clock, before the domain (its climate checks the
    ! indices it names) and the startup (the initial climate carries the same
    ! index values as the time loop). On a restart run, resume them from the
    ! bundle.
    call indices_init(idx, path_par, ts%time_rel)
    if (trim(restart) /= "None") call indices_restart_read(idx, trim(restart))

    ! Initialize the domain: sub-models + hi-res hub + coupler maps. The domain
    ! reads the timeline values it needs from the same [run_step] group.
    call domain_init(dom, path_par, ts%time, idx, timeline_group=trim(run_step))

    ! Define regions of interest for 1D output (must precede the first yelmo_update).
    call domain_regions_init(dom, trim(outfldr))

    ! Cold start: build the initial boundary state. Restart: restore the bundle
    ! (incl. the shared bsl), rebuild the hi-res hub from the restored models,
    ! then re-establish the climate/smb and marine-shelf forcing from the
    ! restored state (the bundle does not hold them), so the first step and the
    ! first output see a valid boundary state.
    call domain_startup(dom, ts, bsl, trim(restart), idx)
    if (trim(restart) /= "None") then
        call step_climate(dom, ts, idx)
        call step_surface(dom, ts)
        call step_shelf(dom, ts)
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
    call domain_write_cmip(dom, trim(outfldr), ts, init=.TRUE.)

    ! === main time loop ===
    ! Output and restarts are written at the top of the loop for the current
    ! time (time_init on the first pass), then the loop exits once the run is
    ! finished (time_end reached, or the kill switch tripped), else the time is
    ! advanced and the domain stepped. The final state and a final restart
    ! bundle are always written.
    call timer_init(tmr, [character(len=16) :: "output", "forcing", "relax", "optimize", &
                          "isostasy", "couple_yelmo", "icesheet", "couple_htopo", &
                          "climate", "surface", "shelf"])

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
            call indices_write_step(idx, dom%yelmo%reg%fnm, ts%time)
            call timer_write_totals(tmr, ts%time, "s", trim(outfldr)//"timing.txt", &
                                    init=timing_init)
            timing_init = .FALSE.
        end if

        ! CMIP-formatted output ([output] write_cmip, every dt_cmip).
        call domain_write_cmip(dom, trim(outfldr), ts)

        ! Restart bundle (domain + shared bsl + indices).
        if (timeout_check(tm_rst, ts%time) .or. ts%is_finished) then
            call run_restart_write(dom, bsl, ts%time, idx)
        end if

        call timer_step(tmr, tt_output)

        if (ts%is_finished) exit

        call tstep_update(ts, dtt)
        call tstep_print(ts)

        ! Shared sea level: update once per step, before the domain advances.
        call bsl_update(bsl, ts%time_rel)

        ! Indices: advance every series once per step on the forcing clock
        ! (feedback methods need the response-derivative window); response
        ! variable = ice volume [Gt].
        fvar = dom%yelmo%reg%V_ice * dom%yelmo%bnd%c%rho_ice * 1e-3_wp
        call indices_update(idx, ts%time_rel, var=fvar)
        call timer_step(tmr, tt_forcing)

        ! === coupling sequence ===
        call step_relax(dom, ts)          ! topography relaxation (relax)
        call timer_step(tmr, tt_relax)
        call step_optimize(dom, ts)       ! cb_ref/tf_corr optimization (opt)
        call timer_step(tmr, tt_optimize)
        call step_isostasy(dom, ts, bsl)  ! bedrock + sea level, this step
        call timer_step(tmr, tt_isostasy)
        call couple_to_yelmo(dom)         ! bedrock now; smb + shelf melt lag one step
        call step_ice_bnd(dom, ts, idx)   ! prescribed ice in the fixed/relaxation cells
        call timer_step(tmr, tt_couple_yelmo)
        call step_icesheet(dom, ts)       ! yelmo_update
        call timer_step(tmr, tt_icesheet)
        call couple_yelmo_to_htopo(dom)   ! hi-res geometry from the models
        call timer_step(tmr, tt_couple_htopo)
        call step_climate(dom, ts, idx)        ! climate (dt_clim cadence)
        call timer_step(tmr, tt_climate)
        call step_surface(dom, ts)        ! surface mass balance + temperature
        call timer_step(tmr, tt_surface)
        call step_shelf(dom, ts)          ! shelf-base melt + temperature
        call timer_step(tmr, tt_shelf)

        ! Forcing-increment restart each time an index moves by more than its
        ! restart_every_df (folders restart-<n>), so a ramp can be branched at
        ! fixed forcing levels.
        if (indices_restart_due(idx)) then
            write(*,*) "yelmox: forcing-increment restart "//trim(indices_restart_fldr(idx))// &
                       " at time =", ts%time
            call run_restart_write(dom, bsl, ts%time, idx, fldr=trim(indices_restart_fldr(idx)))
        end if

        ! Kill switch: stop once the response has equilibrated at a forcing
        ! bound (any index). The top of the loop then writes the final state and exits.
        if (indices_kill(idx)) then
            write(*,*) "yelmox: index kill switch tripped at time =", ts%time
            ts%is_finished = .TRUE.
        end if
    end do

    write(*,*)
    write(*,*) "yelmox: run complete at time =", ts%time
    write(*,*) "  H_ice max   =", maxval(dom%yelmo%tpo%now%H_ice)
    write(*,*)
    write(*,*) "yelmox: wall time of the time loop by component"
    call timer_print_totals(tmr, "s")

    ! Finalize Yelmo (deallocates model state) -- must come after the last
    ! access to dom%yelmo (it deallocates tpo%now%H_ice etc).
    call yelmo_end(dom%yelmo, time=ts%time)
    call climate_end(dom%cl)

end program yelmox
