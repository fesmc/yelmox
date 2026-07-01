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
    use timeout
    use yelmo, only : yelmo_load_command_line_args, wp, &
                      yelmo_write_init, yelmo_write_step, yelmo_end
    use htopo, only : htopo_write_init, htopo_write_step
    use yelmox_domain

    implicit none

    character(len=512) :: path_par
    type(tstep_class)  :: ts
    type(ice_domain)   :: dom
    type(timeout_class) :: tm_2D

    character(len=512) :: file2D, file2D_topo
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

    ! Cold start: build the initial boundary state. Restart: restore the bundle
    ! and rebuild the hi-res hub from the restored models.
    if (trim(dom%ctl%restart) == "None") then
        call domain_init_state(dom, ts)
    else
        call domain_restart_read(dom, trim(dom%ctl%restart), ts%time)
        call refresh_htopo(dom)
    end if

    write(*,*)
    write(*,*) "yelmox_mg: domain initialized"
    write(*,*) "  domain      : "//trim(dom%ctl%domain)
    write(*,*) "  Yelmo grid  : "//trim(dom%ctl%grid_yelmo), dom%yelmo%grd%nx, dom%yelmo%grd%ny
    write(*,*) "  topo grid   : "//trim(dom%ctl%grid_name),  dom%topo%nx,      dom%topo%ny
    write(*,*) "  coupler maps: ", dom%cpl%nmaps
    write(*,*)

    ! === output setup (2D; run from the output folder, yelmox convention) ===
    file2D      = "yelmo2D.nc"
    file2D_topo = "htopo2D.nc"     ! hi-res reference geometry over time
    call timeout_init(tm_2D, path_par, "tm_2D", "heavy", time_init, time_end)
    if (tm_2D%active) then
        call yelmo_write_init(dom%yelmo, file2D, time_init=ts%time, units="years")
        call yelmo_write_step(dom%yelmo, file2D, ts%time, compare_pd=.FALSE.)
        call htopo_write_init(dom%topo, file2D_topo, time_init=ts%time)
        call htopo_write_step(dom%topo, file2D_topo, ts%time)
    end if

    ! === main time loop ===
    call tstep_print_header(ts)
    do while (.not. ts%is_finished)
        call tstep_update(ts, dom%ctl%dtt)
        call tstep_print(ts)

        call yelmox_step(dom, ts)

        if (tm_2D%active .and. timeout_check(tm_2D, ts%time)) then
            call yelmo_write_step(dom%yelmo, file2D, ts%time, compare_pd=.FALSE.)
            call htopo_write_step(dom%topo, file2D_topo, ts%time)
        end if

        if (dom%ctl%dt_restart > 0.0_wp .and. &
            mod(nint(ts%time*100), nint(dom%ctl%dt_restart*100)) == 0) then
            call domain_restart_write(dom, ts%time)
        end if
    end do

    ! Always capture the final state + a final restart bundle.
    if (tm_2D%active) then
        call yelmo_write_step(dom%yelmo, file2D, ts%time, compare_pd=.FALSE.)
        call htopo_write_step(dom%topo, file2D_topo, ts%time)
    end if
    call domain_restart_write(dom, ts%time)

    ! NOTE: yelmo_end is deferred until the region-of-interest setup
    ! (yelmo_regions_init) is lifted in — it finalizes those structures. ncio
    ! closes each write, so yelmo2D.nc is complete without it.

    write(*,*)
    write(*,*) "yelmox_mg: run complete at time =", ts%time
    write(*,*) "  H_ice max   =", maxval(dom%yelmo%tpo%now%H_ice)

end program yelmox_mg
