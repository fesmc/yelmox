program yelmox_mg
    ! Multigrid yelmox driver.
    !
    ! Bring-up stage: initialize one ice_domain (all sub-models on the Yelmo grid
    ! for now) plus the hi-res topography reference hub and the coupler maps,
    ! build the initial boundary state, and run the coupling time loop. No
    ! NetCDF output yet -- output + parity checks against yelmox.f90 come next.
    ! Lifted and adapted from yelmox.f90 (see docs/multigrid.md).

    use nml
    use timestepping
    use yelmo, only : yelmo_load_command_line_args, wp
    use yelmox_domain

    implicit none

    character(len=512) :: path_par
    type(tstep_class)  :: ts
    type(ice_domain)   :: dom

    character(len=56)  :: tstep_method
    real(wp)           :: tstep_const, time_init, time_end

    ! Parameter file path from the command line (runme passes it per run).
    call yelmo_load_command_line_args(path_par)

    ! Timestepping (driver-owned; shared across domains in bipolar runs).
    call nml_read(path_par, "ctrl", "tstep_method", tstep_method)
    call nml_read(path_par, "ctrl", "tstep_const",  tstep_const)
    call nml_read(path_par, "ctrl", "time_init",    time_init)
    call nml_read(path_par, "ctrl", "time_end",     time_end)
    call tstep_init(ts, time_init, time_end, method=tstep_method, units="year", &
                    time_ref=1950.0_wp, const_rel=tstep_const)

    ! Initialize the domain: sub-models + hi-res hub + coupler maps.
    call domain_init(dom, path_par, ts%time, ts%time_rel)

    ! Build the initial boundary state and Yelmo state variables.
    call domain_init_state(dom, ts)

    write(*,*)
    write(*,*) "yelmox_mg: domain initialized"
    write(*,*) "  domain      : "//trim(dom%ctl%domain)
    write(*,*) "  Yelmo grid  : "//trim(dom%ctl%grid_yelmo), dom%yelmo%grd%nx, dom%yelmo%grd%ny
    write(*,*) "  topo grid   : "//trim(dom%ctl%grid_name),  dom%topo%nx,      dom%topo%ny
    write(*,*) "  coupler maps: ", dom%cpl%nmaps
    write(*,*)

    ! === main time loop (no output yet) ===
    call tstep_print_header(ts)
    do while (.not. ts%is_finished)
        call tstep_update(ts, dom%ctl%dtt)
        call tstep_print(ts)
        call yelmox_step(dom, ts)
    end do

    write(*,*)
    write(*,*) "yelmox_mg: run complete at time =", ts%time
    write(*,*) "  H_ice max   =", maxval(dom%yelmo%tpo%now%H_ice)

end program yelmox_mg
