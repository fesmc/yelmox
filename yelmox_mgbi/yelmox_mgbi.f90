program yelmox_mgbi
    ! Multigrid yelmox driver -- bipolar / multi-domain.
    !
    ! Advances several ice domains on a shared timeline. The single command-line
    ! argument is a control file that provides the shared timestepping + output
    ! cadence ([ctrl], [tm_1D], [tm_2D]) and lists the per-domain parameter files
    ! ([domains] par_files). Each listed file is a full single-domain nml (as used
    ! by yelmox_mg); each domain is a self-contained ice_domain (own grids,
    ! coupler and sub-models) and writes to a subfolder named after its domain.
    !
    ! All per-domain physics + coupling live in libs/yelmox_domain.f90, shared
    ! with the single-domain driver (yelmox_mg); this program only adds the
    ! control-file parsing and the loop over domains. See docs/multigrid.md.

    use nml
    use timestepping
    use timeout
    use yelmo, only : yelmo_load_command_line_args, wp, yelmo_end
    use yelmox_domain

    implicit none

    integer, parameter :: NDOM_MAX = 8

    character(len=512) :: path_ctrl
    character(len=512) :: par_files(NDOM_MAX)
    character(len=512), allocatable :: path_par(:), outfldr(:)
    type(tstep_class)               :: ts
    type(ice_domain), allocatable   :: dom(:)
    type(timeout_class)             :: tm_2D, tm_1D
    logical                         :: do_2D, do_1D
    integer                         :: nd, k

    character(len=56) :: tstep_method
    real(wp)          :: tstep_const, time_init, time_end, dtt

    ! Control file from the command line (runme stages the domain files alongside).
    call yelmo_load_command_line_args(path_ctrl)

    ! Shared timestepping (driver-owned) from the control file.
    call nml_read(path_ctrl, "ctrl", "tstep_method", tstep_method)
    call nml_read(path_ctrl, "ctrl", "tstep_const",  tstep_const)
    call nml_read(path_ctrl, "ctrl", "time_init",    time_init)
    call nml_read(path_ctrl, "ctrl", "time_end",     time_end)
    call nml_read(path_ctrl, "ctrl", "dtt",          dtt)
    call tstep_init(ts, time_init, time_end, method=tstep_method, units="year", &
                    time_ref=1950.0_wp, const_rel=tstep_const)

    ! Domain list: [domains] par_files = "fileA.nml" "fileB.nml" ...
    par_files = ""
    call nml_read(path_ctrl, "domains", "par_files", par_files)
    nd = count(len_trim(par_files) > 0)
    if (nd < 1) then
        write(*,*) "yelmox_mgbi: error -- [domains] par_files lists no domains."
        stop 1
    end if
    allocate(path_par(nd), outfldr(nd), dom(nd))
    path_par(1:nd) = par_files(1:nd)

    ! === Per-domain initialization ===
    do k = 1, nd

        call domain_init(dom(k), trim(path_par(k)), ts%time, ts%time_rel)

        ! Each domain writes to a subfolder named after its domain.
        outfldr(k) = trim(dom(k)%ctl%domain)//"/"
        call execute_command_line('mkdir -p "'//trim(outfldr(k))//'"')

        call domain_regions_init(dom(k), trim(outfldr(k)))

        if (trim(dom(k)%ctl%restart) == "None") then
            call domain_init_state(dom(k), ts)
        else
            call domain_restart_read(dom(k), trim(dom(k)%ctl%restart), ts)
            call refresh_htopo(dom(k))
        end if

        write(*,*)
        write(*,*) "yelmox_mgbi: domain initialized"
        write(*,*) "  domain      : "//trim(dom(k)%ctl%domain)
        write(*,*) "  Yelmo grid  : "//trim(dom(k)%ctl%grid_yelmo), dom(k)%yelmo%grd%nx, dom(k)%yelmo%grd%ny
        write(*,*) "  topo grid   : "//trim(dom(k)%ctl%grid_name),  dom(k)%topo%nx,      dom(k)%topo%ny
        write(*,*) "  coupler maps: ", dom(k)%cpl%nmaps
        write(*,*) "  output dir  : "//trim(outfldr(k))
        write(*,*)

    end do

    ! === Output setup (shared cadence, from the control file) ===
    call timeout_init(tm_2D, path_ctrl, "tm_2D", "heavy", time_init, time_end)
    call timeout_init(tm_1D, path_ctrl, "tm_1D", "small", time_init, time_end)
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
    write(*,*) "yelmox_mgbi: run complete at time =", ts%time
    do k = 1, nd
        write(*,*) "  "//trim(dom(k)%ctl%domain)//" H_ice max =", maxval(dom(k)%yelmo%tpo%now%H_ice)
    end do

    do k = 1, nd
        call yelmo_end(dom(k)%yelmo, time=ts%time)
    end do

end program yelmox_mgbi
