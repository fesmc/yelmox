program yelmox_mg
    ! Multigrid yelmox driver (single domain).
    !
    ! Initializes one ice_domain (each sub-model on its own configurable grid)
    ! plus the hi-res topography reference hub and the coupler maps, builds the
    ! initial boundary state (or restores a restart bundle), and runs the coupling
    ! time loop with per-module output. The multi-domain (bipolar) variant lives
    ! in yelmox_mgbi/. See docs/multigrid.md and libs/yelmox_domain.f90.

    use nml
    use timestepping
    use timeout
    use yelmo, only : yelmo_load_command_line_args, wp, yelmo_end
    use yelmox_domain

    implicit none

    character(len=512) :: path_par
    type(tstep_class)  :: ts
    type(ice_domain)   :: dom
    type(timeout_class) :: tm_2D, tm_1D

    character(len=512) :: outfldr
    character(len=56)  :: tstep_method
    real(wp)           :: tstep_const, time_init, time_end, dtt

    ! Parameter file path from the command line (runme passes it per run).
    call yelmo_load_command_line_args(path_par)

    ! Timestepping (driver-owned; the [ctrl] group holds the shared timeline).
    call nml_read(path_par, "ctrl", "tstep_method", tstep_method)
    call nml_read(path_par, "ctrl", "tstep_const",  tstep_const)
    call nml_read(path_par, "ctrl", "time_init",    time_init)
    call nml_read(path_par, "ctrl", "time_end",     time_end)
    call nml_read(path_par, "ctrl", "dtt",          dtt)
    call tstep_init(ts, time_init, time_end, method=tstep_method, units="year", &
                    time_ref=1950.0_wp, const_rel=tstep_const)

    ! Single-domain runs write to the run dir.
    outfldr = "./"

    ! Initialize the domain: sub-models + hi-res hub + coupler maps.
    call domain_init(dom, path_par, ts%time, ts%time_rel)

    ! Inject the driver-owned timeline values the domain logic needs.
    dom%ctl%tstep_method = tstep_method
    dom%ctl%dtt          = dtt

    ! Define regions of interest for 1D output (must precede the first yelmo_update).
    call domain_regions_init(dom, trim(outfldr))

    ! Cold start: build the initial boundary state. Restart: restore the bundle
    ! and rebuild the hi-res hub from the restored models.
    if (trim(dom%ctl%restart) == "None") then
        call domain_init_state(dom, ts)
    else
        call domain_restart_read(dom, trim(dom%ctl%restart), ts)
        call refresh_htopo(dom)
    end if

    write(*,*)
    write(*,*) "yelmox_mg: domain initialized"
    write(*,*) "  domain      : "//trim(dom%ctl%domain)
    write(*,*) "  Yelmo grid  : "//trim(dom%ctl%grid_yelmo), dom%yelmo%grd%nx, dom%yelmo%grd%ny
    write(*,*) "  topo grid   : "//trim(dom%ctl%grid_name),  dom%topo%nx,      dom%topo%ny
    write(*,*) "  coupler maps: ", dom%cpl%nmaps
    write(*,*)

    ! === output setup (2D; one file per module, on its own grid) ===
    call timeout_init(tm_2D, path_par, "tm_2D", "heavy", time_init, time_end)
    if (tm_2D%active) then
        call domain_write_init(dom, trim(outfldr), ts%time)
        call domain_write_step(dom, trim(outfldr), ts%time)
    end if

    ! === output setup (1D timeseries) ===
    call timeout_init(tm_1D, path_par, "tm_1D", "small", time_init, time_end)
    if (tm_1D%active) then
        call domain_write_1D(dom, trim(outfldr), ts%time, init=.TRUE.)
    end if

    ! === main time loop ===
    call tstep_print_header(ts)
    do while (.not. ts%is_finished)
        call tstep_update(ts, dom%ctl%dtt)
        call tstep_print(ts)

        call yelmox_step(dom, ts)

        if (tm_2D%active .and. timeout_check(tm_2D, ts%time)) then
            call domain_write_step(dom, trim(outfldr), ts%time)
        end if

        if (tm_1D%active .and. timeout_check(tm_1D, ts%time)) then
            call domain_write_1D(dom, trim(outfldr), ts%time)
        end if

        if (dom%ctl%dt_restart > 0.0_wp .and. &
            mod(nint(ts%time*100), nint(dom%ctl%dt_restart*100)) == 0) then
            call domain_restart_write(dom, ts%time)
        end if
    end do

    ! Always capture the final state + a final restart bundle.
    if (tm_2D%active) call domain_write_step(dom, trim(outfldr), ts%time)
    if (tm_1D%active) call domain_write_1D(dom, trim(outfldr), ts%time)
    call domain_restart_write(dom, ts%time)

    write(*,*)
    write(*,*) "yelmox_mg: run complete at time =", ts%time
    write(*,*) "  H_ice max   =", maxval(dom%yelmo%tpo%now%H_ice)

    ! Finalize Yelmo (deallocates model state) -- must come after the last
    ! access to dom%yelmo (it deallocates tpo%now%H_ice etc).
    call yelmo_end(dom%yelmo, time=ts%time)

end program yelmox_mg
