program yelmox_mg
    ! Multigrid yelmox driver.
    !
    ! Advances one or more ice domains on a shared timeline. Each domain is a
    ! self-contained ice_domain (its own grids, coupler and sub-models) built from
    ! its own parameter file: one command-line argument per domain. A single par
    ! file runs one domain (writing to the run dir); several par files run a
    ! multi-domain (e.g. bipolar) setup, each domain writing to a subfolder named
    ! after its domain. See docs/multigrid.md and libs/yelmox_domain.f90.

    use nml
    use timestepping
    use timeout
    use yelmo, only : wp, yelmo_end
    use yelmox_domain

    implicit none

    integer :: nd, k
    character(len=512), allocatable :: path_par(:), outfldr(:)
    type(tstep_class)               :: ts
    type(ice_domain), allocatable   :: dom(:)
    type(timeout_class)             :: tm_2D, tm_1D
    logical                         :: do_2D, do_1D

    character(len=56) :: tstep_method
    real(wp)          :: tstep_const, time_init, time_end

    ! Domains = one per command-line parameter file (>= 1; more than one is a
    ! multi-domain / bipolar run).
    nd = command_argument_count()
    if (nd < 1) then
        write(*,*) "yelmox_mg: error -- provide at least one parameter file."
        stop 1
    end if
    allocate(path_par(nd), outfldr(nd), dom(nd))
    do k = 1, nd
        call get_command_argument(k, path_par(k))
    end do

    ! Timestepping (driver-owned; shared across domains). Read from the first par.
    call nml_read(path_par(1), "ctrl", "tstep_method", tstep_method)
    call nml_read(path_par(1), "ctrl", "tstep_const",  tstep_const)
    call nml_read(path_par(1), "ctrl", "time_init",    time_init)
    call nml_read(path_par(1), "ctrl", "time_end",     time_end)
    call tstep_init(ts, time_init, time_end, method=tstep_method, units="year", &
                    time_ref=1950.0_wp, const_rel=tstep_const)

    ! === Per-domain initialization ===
    do k = 1, nd

        ! Initialize the domain: sub-models + hi-res hub + coupler maps.
        call domain_init(dom(k), path_par(k), ts%time, ts%time_rel)

        ! Output folder: run dir for a single domain, per-domain subfolder otherwise.
        if (nd == 1) then
            outfldr(k) = "./"
        else
            outfldr(k) = trim(dom(k)%ctl%domain)//"/"
            call execute_command_line('mkdir -p "'//trim(outfldr(k))//'"')
        end if

        ! Regions of interest for 1D output (must precede the first yelmo_update).
        call domain_regions_init(dom(k), trim(outfldr(k)))

        ! Cold start: build the initial boundary state. Restart: restore the bundle
        ! and rebuild the hi-res hub from the restored models.
        if (trim(dom(k)%ctl%restart) == "None") then
            call domain_init_state(dom(k), ts)
        else
            call domain_restart_read(dom(k), trim(dom(k)%ctl%restart), ts)
            call refresh_htopo(dom(k))
        end if

        write(*,*)
        write(*,*) "yelmox_mg: domain initialized"
        write(*,*) "  domain      : "//trim(dom(k)%ctl%domain)
        write(*,*) "  Yelmo grid  : "//trim(dom(k)%ctl%grid_yelmo), dom(k)%yelmo%grd%nx, dom(k)%yelmo%grd%ny
        write(*,*) "  topo grid   : "//trim(dom(k)%ctl%grid_name),  dom(k)%topo%nx,      dom(k)%topo%ny
        write(*,*) "  coupler maps: ", dom(k)%cpl%nmaps
        write(*,*) "  output dir  : "//trim(outfldr(k))
        write(*,*)

    end do

    ! === Output setup (shared cadence, from the first par) ===
    call timeout_init(tm_2D, path_par(1), "tm_2D", "heavy", time_init, time_end)
    call timeout_init(tm_1D, path_par(1), "tm_1D", "small", time_init, time_end)
    do k = 1, nd
        if (tm_2D%active) then
            call domain_write_init(dom(k), trim(outfldr(k)), ts%time)
            call domain_write_step(dom(k), trim(outfldr(k)), ts%time)
        end if
        if (tm_1D%active) call domain_write_1D(dom(k), trim(outfldr(k)), ts%time, init=.TRUE.)
    end do

    ! === Main time loop ===
    call tstep_print_header(ts)
    do while (.not. ts%is_finished)
        call tstep_update(ts, dom(1)%ctl%dtt)
        call tstep_print(ts)

        ! Advance every domain by one coupling step.
        do k = 1, nd
            call yelmox_step(dom(k), ts)
        end do

        ! Evaluate the shared output cadence once (timeout_check advances state).
        do_2D = tm_2D%active .and. timeout_check(tm_2D, ts%time)
        do_1D = tm_1D%active .and. timeout_check(tm_1D, ts%time)

        do k = 1, nd
            if (do_2D) call domain_write_step(dom(k), trim(outfldr(k)), ts%time)
            if (do_1D) call domain_write_1D(dom(k), trim(outfldr(k)), ts%time)
            if (dom(k)%ctl%dt_restart > 0.0_wp .and. &
                mod(nint(ts%time*100), nint(dom(k)%ctl%dt_restart*100)) == 0) then
                call domain_restart_write(dom(k), ts%time, outfldr=trim(outfldr(k)))
            end if
        end do
    end do

    ! === Finalize: capture the final state + a final restart bundle per domain ===
    do k = 1, nd
        if (tm_2D%active) call domain_write_step(dom(k), trim(outfldr(k)), ts%time)
        if (tm_1D%active) call domain_write_1D(dom(k), trim(outfldr(k)), ts%time)
        call domain_restart_write(dom(k), ts%time, outfldr=trim(outfldr(k)))
    end do

    write(*,*)
    write(*,*) "yelmox_mg: run complete at time =", ts%time
    do k = 1, nd
        write(*,*) "  "//trim(dom(k)%ctl%domain)//" H_ice max =", maxval(dom(k)%yelmo%tpo%now%H_ice)
    end do

    ! Finalize Yelmo (deallocates model state) -- after the last access to dom%yelmo.
    do k = 1, nd
        call yelmo_end(dom(k)%yelmo, time=ts%time)
    end do

end program yelmox_mg
