program yelmox
    ! Multigrid yelmox driver (single domain).
    !
    ! Initializes one ice_domain (each sub-model on its own configurable grid)
    ! plus the hi-res topography reference hub and the coupler maps, builds the
    ! initial boundary state (or restores a restart bundle), and runs the coupling
    ! time loop with per-module output. The multi-domain (bipolar) variant lives
    ! in yelmox_bipolar/. See docs/multigrid.md and libs/yelmox_domain.f90.

    use nml
    use timestepping
    use timeout
    use yelmo, only : yelmo_load_command_line_args, wp, yelmo_end
    use fastisostasy, only : bsl_class, bsl_init, bsl_update
    use tsgen, only : tsgen_class, tsgen_init, tsgen_update
    use yelmox_domain

    implicit none

    character(len=512) :: path_par
    type(tstep_class)  :: ts
    type(ice_domain)   :: dom
    type(bsl_class)    :: bsl        ! shared, driver-owned barystatic sea level
    type(timeout_class) :: tm_2D, tm_1D

    character(len=512) :: outfldr
    real(wp)           :: dtt

    ! Transient time-series forcing (tsgen), owned by the driver. The single
    ! forcing value f_now is mapped onto the snapclim anomalies via per-channel
    ! gains ([tsforcing]): dTa = f_now*f_ta, dTo = f_now*f_to, dSo = f_now*f_so.
    type(tsgen_class) :: tsg
    logical  :: tsf_active
    real(wp) :: tsf_f_ta, tsf_f_to, tsf_f_so
    real(wp) :: fvar

    ! Parameter file path from the command line (runme passes it per run).
    call yelmo_load_command_line_args(path_par)

    ! Timestepping (driver-owned; the [ctrl] group holds the shared timeline).
    call timeline_init(ts, dtt, path_par, "ctrl")

    ! Single-domain runs write to the run dir.
    outfldr = "./"

    ! Shared, driver-owned barystatic sea level (one per run).
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, ts%time_rel)

    ! Initialize the domain: sub-models + hi-res hub + coupler maps. The domain
    ! reads the timeline values it needs from the same [ctrl] group.
    call domain_init(dom, path_par, ts%time)

    ! Define regions of interest for 1D output (must precede the first yelmo_update).
    call domain_regions_init(dom, trim(outfldr))

    ! Transient time-series forcing (tsgen -> snapclim anomalies). Initialize
    ! before startup so the initial (cold-start) climate carries the same
    ! anomalies as the time loop. tsgen reads its own parameters from [tsgen];
    ! tsgen_init sets f_now to the series value at the initial time.
    call nml_read(path_par, "tsforcing", "active", tsf_active)
    call nml_read(path_par, "tsforcing", "f_ta",   tsf_f_ta)
    call nml_read(path_par, "tsforcing", "f_to",   tsf_f_to)
    call nml_read(path_par, "tsforcing", "f_so",   tsf_f_so)
    if (tsf_active) then
        call tsgen_init(tsg, path_par, ts%time)
        write(*,*) "yelmox: transient forcing active (tsgen), f_ta/f_to/f_so =", &
                    tsf_f_ta, tsf_f_to, tsf_f_so
    end if

    ! Cold start: build the initial boundary state. Restart: restore the bundle
    ! (incl. the shared bsl) and rebuild the hi-res hub from the restored models.
    if (tsf_active) then
        call domain_startup(dom, ts, bsl, dTa=tsg%f_now*tsf_f_ta, &
                            dTo=tsg%f_now*tsf_f_to, dSo=tsg%f_now*tsf_f_so)
    else
        call domain_startup(dom, ts, bsl)
    end if

    write(*,*)
    write(*,*) "yelmox: domain initialized"
    write(*,*) "  domain      : "//trim(dom%ctl%domain)
    write(*,*) "  Yelmo grid  : "//trim(dom%ctl%grid_yelmo), dom%yelmo%grd%G%nx, dom%yelmo%grd%G%ny
    write(*,*) "  topo grid   : "//trim(dom%ctl%grid_name),  dom%topo%nx,      dom%topo%ny
    write(*,*) "  coupler maps: ", dom%cpl%nmaps
    write(*,*)

    ! === output setup (2D; one file per module, on its own grid) ===
    call timeout_init(tm_2D, path_par, "tm_2D", "heavy", ts%time_init, ts%time_end)
    if (tm_2D%active) then
        call domain_write_init(dom, trim(outfldr), ts%time)
        call domain_write_step(dom, trim(outfldr), ts%time)
    end if

    ! === output setup (1D timeseries) ===
    call timeout_init(tm_1D, path_par, "tm_1D", "small", ts%time_init, ts%time_end)
    if (tm_1D%active) then
        call domain_write_1D(dom, trim(outfldr), ts%time, init=.TRUE.)
    end if

    ! === main time loop ===
    call tstep_print_header(ts)
    do while (.not. ts%is_finished)
        call tstep_update(ts, dtt)
        call tstep_print(ts)

        ! Shared sea level: update once per step, before the domain advances.
        call bsl_update(bsl, ts%time_rel)

        if (tsf_active) then
            ! Advance the forcing series every step (feedback methods need the
            ! response-derivative window); response variable = ice volume [Gt].
            fvar = dom%yelmo%reg%V_ice * dom%yelmo%bnd%c%rho_ice * 1e-3_wp
            call tsgen_update(tsg, ts%time, var=fvar)
            call yelmox_step(dom, ts, bsl, dTa=tsg%f_now*tsf_f_ta, &
                             dTo=tsg%f_now*tsf_f_to, dSo=tsg%f_now*tsf_f_so)
        else
            call yelmox_step(dom, ts, bsl)
        end if

        if (tm_2D%active .and. timeout_check(tm_2D, ts%time)) then
            call domain_write_step(dom, trim(outfldr), ts%time)
        end if

        if (tm_1D%active .and. timeout_check(tm_1D, ts%time)) then
            call domain_write_1D(dom, trim(outfldr), ts%time)
        end if

        if (tstep_due(ts%time, dom%ctl%dt_restart)) then
            call run_restart_write(dom, bsl, ts%time)
        end if
    end do

    ! Always capture the final state + a final restart bundle (incl. shared bsl).
    if (tm_2D%active) call domain_write_step(dom, trim(outfldr), ts%time)
    if (tm_1D%active) call domain_write_1D(dom, trim(outfldr), ts%time)
    call run_restart_write(dom, bsl, ts%time)

    write(*,*)
    write(*,*) "yelmox: run complete at time =", ts%time
    write(*,*) "  H_ice max   =", maxval(dom%yelmo%tpo%now%H_ice)

    ! Finalize Yelmo (deallocates model state) -- must come after the last
    ! access to dom%yelmo (it deallocates tpo%now%H_ice etc).
    call yelmo_end(dom%yelmo, time=ts%time)

end program yelmox
