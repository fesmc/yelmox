program yelmox_mgbi
    ! Multigrid yelmox driver -- bipolar / multi-domain.
    !
    ! Advances several ice domains on a shared timeline. The single command-line
    ! argument is one parameter file holding every domain, following the original
    ! yelmox_bipolar convention: each domain's groups carry a domain suffix
    ! (yelmo_north, coupling_south, snap_north, ...), while shared blocks ([ctrl],
    ! [barysealevel], the yelmo physics groups ydyn/ytopo/...) have no suffix.
    ! Keeping every group name distinct also lets `runme -p group.name=val` target
    ! a single domain unambiguously.
    !
    ! [domains] names lists the per-domain suffix tags; each domain is a
    ! self-contained ice_domain (own grids, coupler and sub-models) written to a
    ! subfolder named after its domain. All per-domain physics + coupling live in
    ! libs/yelmox_domain.f90, shared with the single-domain driver (yelmox_mg);
    ! this program only adds the domain loop. See docs/multigrid.md.

    use nml
    use timestepping
    use timeout
    use yelmo, only : yelmo_load_command_line_args, wp, yelmo_end
    use fastisostasy, only : bsl_class, bsl_init, bsl_update, &
                             bsl_restart_read, bsl_restart_write
    use yelmox_domain

    implicit none

    integer, parameter :: NDOM_MAX = 8

    character(len=512) :: path_par
    character(len=64)  :: names(NDOM_MAX)
    character(len=512), allocatable :: outfldr(:)
    type(tstep_class)               :: ts
    type(ice_domain), allocatable   :: dom(:)
    type(bsl_class)                 :: bsl        ! shared, driver-owned sea level
    character(len=512)              :: restart_bsl
    type(timeout_class)             :: tm_2D, tm_1D
    logical                         :: do_2D, do_1D, do_restart
    integer                         :: nd, k

    character(len=56) :: tstep_method
    real(wp)          :: tstep_const, time_init, time_end, dtt

    ! Single parameter file from the command line (holds all domains).
    call yelmo_load_command_line_args(path_par)

    ! Shared timestepping (driver-owned) from the [ctrl] group.
    call nml_read(path_par, "ctrl", "tstep_method", tstep_method)
    call nml_read(path_par, "ctrl", "tstep_const",  tstep_const)
    call nml_read(path_par, "ctrl", "time_init",    time_init)
    call nml_read(path_par, "ctrl", "time_end",     time_end)
    call nml_read(path_par, "ctrl", "dtt",          dtt)
    call tstep_init(ts, time_init, time_end, method=tstep_method, units="year", &
                    time_ref=1950.0_wp, const_rel=tstep_const)

    ! Domain suffix tags: [domains] names = "north" "south" ... -> group suffix
    ! "_north", "_south", etc.
    names = ""
    call nml_read(path_par, "domains", "names", names)
    nd = count(len_trim(names) > 0)
    if (nd < 1) then
        write(*,*) "yelmox_mgbi: error -- [domains] names lists no domains."
        stop 1
    end if
    allocate(outfldr(nd), dom(nd))

    ! === Shared, driver-owned barystatic sea level (one per run, common to every
    !     domain, exactly as in yelmox_bipolar). Restored once from a run-level
    !     bsl_restart.nc ([ctrl] restart_bsl) when restarting. ===
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, ts%time_rel)
    restart_bsl = "None"
    call nml_read(path_par, "ctrl", "restart_bsl", restart_bsl)
    if (trim(restart_bsl) /= "None") then
        call bsl_restart_read(bsl, trim(restart_bsl)//"/bsl_restart.nc")
        call bsl_update(bsl, ts%time_rel)
    end if

    ! === Per-domain initialization ===
    do k = 1, nd

        call domain_init(dom(k), path_par, ts%time, &
                         group_suffix="_"//trim(names(k)))

        ! Inject the driver-owned timeline values the domain logic needs.
        dom(k)%ctl%tstep_method = tstep_method
        dom(k)%ctl%dtt          = dtt

        ! Each domain writes to a subfolder named after its domain.
        outfldr(k) = trim(dom(k)%ctl%domain)//"/"
        call execute_command_line('mkdir -p "'//trim(outfldr(k))//'"')

        call domain_regions_init(dom(k), trim(outfldr(k)))

        if (trim(dom(k)%ctl%restart) == "None") then
            call domain_init_state(dom(k), ts, bsl)
        else
            call domain_restart_read(dom(k), trim(dom(k)%ctl%restart), ts, bsl)
            call refresh_htopo(dom(k))
        end if

        write(*,*)
        write(*,*) "yelmox_mgbi: domain initialized ("//trim(names(k))//")"
        write(*,*) "  domain      : "//trim(dom(k)%ctl%domain)
        write(*,*) "  Yelmo grid  : "//trim(dom(k)%ctl%grid_yelmo), dom(k)%yelmo%grd%nx, dom(k)%yelmo%grd%ny
        write(*,*) "  topo grid   : "//trim(dom(k)%ctl%grid_name),  dom(k)%topo%nx,      dom(k)%topo%ny
        write(*,*) "  coupler maps: ", dom(k)%cpl%nmaps
        write(*,*) "  output dir  : "//trim(outfldr(k))
        write(*,*)

    end do

    ! === Output setup (shared cadence, from [tm_1D]/[tm_2D]) ===
    call timeout_init(tm_2D, path_par, "tm_2D", "heavy", time_init, time_end)
    call timeout_init(tm_1D, path_par, "tm_1D", "small", time_init, time_end)
    do k = 1, nd
        if (tm_2D%active) then
            call domain_write_init(dom(k), trim(outfldr(k)), ts%time)
            call domain_write_step(dom(k), trim(outfldr(k)), ts%time)
        end if
        if (tm_1D%active) call domain_write_1D(dom(k), trim(outfldr(k)), ts%time, init=.TRUE.)
    end do

    ! === Main time loop (shared timeline + dtt) ===
    call tstep_print_header(ts)
    do while (.not. ts%is_finished)
        call tstep_update(ts, dtt)
        call tstep_print(ts)

        ! Shared sea level: update once per step, before any domain advances, so
        ! the two domains' isostasy see (and, under fastiso, contribute to) the
        ! same barystatic sea level.
        call bsl_update(bsl, ts%time_rel)

        do k = 1, nd
            call yelmox_step(dom(k), ts, bsl)
        end do

        ! Evaluate the shared output cadence once (timeout_check advances state).
        do_2D = tm_2D%active .and. timeout_check(tm_2D, ts%time)
        do_1D = tm_1D%active .and. timeout_check(tm_1D, ts%time)

        do_restart = .false.
        do k = 1, nd
            if (do_2D) call domain_write_step(dom(k), trim(outfldr(k)), ts%time)
            if (do_1D) call domain_write_1D(dom(k), trim(outfldr(k)), ts%time)
            if (dom(k)%ctl%dt_restart > 0.0_wp .and. &
                mod(nint(ts%time*100), nint(dom(k)%ctl%dt_restart*100)) == 0) then
                call domain_restart_write(dom(k), ts%time, outfldr=trim(outfldr(k)))
                do_restart = .true.
            end if
        end do

        ! One shared bsl restart at the run root when any domain wrote a bundle.
        if (do_restart) &
            call bsl_restart_write(bsl, trim(restart_bundle_dir(ts%time))//"/bsl_restart.nc", ts%time)
    end do

    ! === Finalize: capture the final state + a final restart bundle per domain ===
    do k = 1, nd
        if (tm_2D%active) call domain_write_step(dom(k), trim(outfldr(k)), ts%time)
        if (tm_1D%active) call domain_write_1D(dom(k), trim(outfldr(k)), ts%time)
        call domain_restart_write(dom(k), ts%time, outfldr=trim(outfldr(k)))
    end do
    call bsl_restart_write(bsl, trim(restart_bundle_dir(ts%time))//"/bsl_restart.nc", ts%time)

    write(*,*)
    write(*,*) "yelmox_mgbi: run complete at time =", ts%time
    do k = 1, nd
        write(*,*) "  "//trim(dom(k)%ctl%domain)//" H_ice max =", maxval(dom(k)%yelmo%tpo%now%H_ice)
    end do

    do k = 1, nd
        call yelmo_end(dom(k)%yelmo, time=ts%time)
    end do

end program yelmox_mgbi
