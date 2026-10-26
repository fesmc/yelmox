program yelmox
    ! Multigrid yelmox driver (single domain).
    !
    ! Initializes one kryos_domain (each sub-model on its own configurable grid)
    ! plus the hi-res topography reference hub and the coupler maps, builds the
    ! initial boundary state (or restores a restart bundle), and runs the coupling
    ! time loop with per-module output. The multi-domain (bipolar) variant lives
    ! in yelmox_bipolar/. See docs/multigrid.md and libs/kryos*.f90.

    use nml
    use ncio,      only : nc_set_lossless
    use timestepping
    use timeout
    use precision, only : dp
    use yelmo, only : yelmo_load_command_line_args, wp, yelmo_end
    use fastisostasy, only : bsl_class, bsl_init, bsl_update
    use kryos,          only : kryos_domain, domain_init, domain_coefs_used
    use yelmox_climate, only : climate_end
    use kryos_regions,  only : domain_regions_init
    use kryos_coupling, only : step_relax, step_optimize, step_isostasy, couple_to_yelmo, &
                               step_ice_bnd, step_icesheet, couple_yelmo_to_htopo, step_climate, &
                               step_surface, step_shelf
    use kryos_startup,  only : domain_startup, run_restart_write
    use kryos_forcing,  only : coefs_class, coefs_init, coefs_update, coefs_kill, &
                               coefs_restart_due, coefs_restart_fldr, &
                               coefs_restart_read, coefs_check_used
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

    ! Write all output lossless ([ctrl]), ignoring the precision of each field.
    logical            :: write_lossless

    ! Coefficients ([forcing] coefs, one group [coef_<name>] each: tsgen, or
    ! derived and set by the model), owned by the driver and advanced once per
    ! step on the forcing time time_rel; read by name by the climate (the
    ! terms) and the prescribed ice. The registry also owns the kill switch
    ! and the forcing-increment restarts (see libs/kryos_forcing.f90).
    type(coefs_class)   :: coefs
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
    ! Optional: write all output lossless (ignore prec_abs/prec_rel), e.g. for bit-for-bit comparisons
    write_lossless = .FALSE.
    if (nml_has_param(path_par, "ctrl", "write_lossless")) &
        call nml_read(path_par, "ctrl", "write_lossless", write_lossless)
    call nc_set_lossless(write_lossless)
    call tstep_init(ts, path_par, trim(run_step), dtt, time_ref=calendar_ref, cal=calendar)

    ! Single-domain runs write to the run dir.
    outfldr = "./"

    ! Shared, driver-owned barystatic sea level (one per run).
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, ts%time_rel)

    ! Coefficients on the forcing time, before the domain (its climate checks
    ! the coefficients it names) and the startup (the initial climate carries
    ! the same coefficient values as the time loop). On a restart run, resume them from the
    ! bundle.
    call coefs_init(coefs, path_par, ts%time_rel)
    if (trim(restart) /= "None") call coefs_restart_read(coefs, trim(restart))

    ! Initialize the domain: sub-models + hi-res hub + coupler maps, at the
    ! start of the timeline ([run_step] group).
    call domain_init(dom, path_par, ts, coefs, timeline_group=trim(run_step))

    ! Every coefficient of [forcing] coefs is read by the domain.
    call coefs_check_used(coefs, domain_coefs_used(dom))

    ! Define regions of interest for 1D output (must precede the first yelmo_update).
    call domain_regions_init(dom, trim(outfldr))

    ! Cold start: build the initial boundary state. Restart: restore the bundle
    ! (incl. the shared bsl) and continue its state: the hi-res hub rebuilt
    ! from the restored models, the climate (and a stateless surface)
    ! re-derived from it (see domain_startup).
    call domain_startup(dom, ts, dtt, bsl, trim(restart), coefs)

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
    if (tm_1D%active)   call domain_write_1D(dom, trim(outfldr), ts%time, coefs, init=.TRUE.)
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
            call domain_write_1D(dom, trim(outfldr), ts%time, coefs)
            call timer_write_totals(tmr, ts%time, "s", trim(outfldr)//"timing.txt", &
                                    init=timing_init)
            timing_init = .FALSE.
        end if

        ! CMIP-formatted output ([output] write_cmip, every dt_cmip).
        call domain_write_cmip(dom, trim(outfldr), ts)

        ! Restart bundle (domain + shared bsl + coefficients).
        if (timeout_check(tm_rst, ts%time) .or. ts%is_finished) then
            call run_restart_write(dom, bsl, ts%time, coefs)
        end if

        call timer_step(tmr, tt_output)

        if (ts%is_finished) exit

        call tstep_update(ts, dtt)
        call tstep_print(ts)

        ! Shared sea level: update once per step, before the domain advances.
        call bsl_update(bsl, ts%time_rel)

        ! Coefficients: advance every tsgen coefficient once per step on the forcing time
        ! (feedback methods need the response-derivative window); response
        ! variable = ice volume [Gt].
        fvar = dom%yelmo%reg%V_ice * dom%yelmo%bnd%c%rho_ice * 1e-3_wp
        call coefs_update(coefs, ts%time_rel, var=fvar)
        call timer_step(tmr, tt_forcing)

        ! === coupling sequence ===
        call step_relax(dom, ts)          ! topography relaxation (relax)
        call timer_step(tmr, tt_relax)
        call step_optimize(dom, ts, dtt)  ! cb_ref/tf_corr optimization (opt)
        call timer_step(tmr, tt_optimize)
        call step_isostasy(dom, ts, bsl)  ! bedrock + sea level, this step
        call timer_step(tmr, tt_isostasy)
        call couple_to_yelmo(dom)         ! bedrock now; smb + shelf melt lag one step
        call step_ice_bnd(dom, ts, coefs)   ! prescribed ice in the fixed/relaxation cells
        call timer_step(tmr, tt_couple_yelmo)
        call step_icesheet(dom, ts)       ! yelmo_update
        call timer_step(tmr, tt_icesheet)
        call couple_yelmo_to_htopo(dom)   ! hi-res geometry from the models
        call timer_step(tmr, tt_couple_htopo)
        call step_climate(dom, ts, dtt, coefs)   ! climate (dt_clim cadence)
        call timer_step(tmr, tt_climate)
        call step_surface(dom, ts)        ! surface mass balance + temperature
        call timer_step(tmr, tt_surface)
        call step_shelf(dom, ts)          ! shelf-base melt + temperature
        call timer_step(tmr, tt_shelf)

        ! Forcing-increment restart each time a coefficient moves by more than its
        ! restart_every_df (folders restart-<n>), so a ramp can be branched at
        ! fixed forcing levels.
        if (coefs_restart_due(coefs)) then
            write(*,*) "yelmox: forcing-increment restart "//trim(coefs_restart_fldr(coefs))// &
                       " at time =", ts%time
            call run_restart_write(dom, bsl, ts%time, coefs, fldr=trim(coefs_restart_fldr(coefs)))
        end if

        ! Kill switch: stop once the response has equilibrated at a forcing
        ! bound (any coefficient). The top of the loop then writes the final state and exits.
        if (coefs_kill(coefs)) then
            write(*,*) "yelmox: coefficient kill switch tripped at time =", ts%time
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
