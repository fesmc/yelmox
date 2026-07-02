program yelmox_mgbi
    ! Multigrid yelmox driver -- bipolar (two hemispheric domains).
    !
    ! Advances a Northern- and a Southern-Hemisphere ice_domain on a shared
    ! timeline, following the original yelmox_bipolar convention: one parameter
    ! file holds both domains, each domain's groups carry a hemisphere suffix
    ! (yelmo_north, coupling_south, snap_north, ...), while shared blocks ([ctrl],
    ! [barysealevel], the yelmo physics groups ydyn/ytopo/...) have no suffix.
    ! Distinct group names also let `runme -p group.name=val` target one
    ! hemisphere unambiguously.
    !
    ! [ctrl] active_north / active_south select which hemispheres run; a bipolar
    ! run is never more than these two domains, so they are held as two explicit
    ! ice_domain variables (not an array) -- the inter-domain ocean coupling is
    ! inherently asymmetric (north <-> obm%fn/thetan/tn, south <-> obm%fs/...).
    ! Each domain is a self-contained ice_domain (own grids, coupler, sub-models)
    ! written to a subfolder named after its domain. All per-domain physics live
    ! in libs/yelmox_domain.f90, shared with the single-domain driver (yelmox_mg).
    ! See docs/multigrid.md.

    use nml
    use timestepping
    use timeout
    use yelmo, only : yelmo_load_command_line_args, wp, yelmo_end
    use fastisostasy, only : bsl_class, bsl_init, bsl_update, &
                             bsl_restart_read, bsl_restart_write
    use yelmox_domain

    implicit none

    character(len=512) :: path_par
    type(tstep_class)  :: ts
    type(bsl_class)    :: bsl                 ! shared, driver-owned sea level
    character(len=512) :: restart_bsl

    ! Two explicit hemispheric domains (bipolar: never more than north + south).
    type(ice_domain)   :: dom_north, dom_south
    character(len=512) :: outfldr_north, outfldr_south
    logical            :: active_north, active_south

    type(timeout_class) :: tm_2D, tm_1D
    logical             :: do_2D, do_1D, do_restart

    character(len=56) :: tstep_method
    real(wp)          :: tstep_const, time_init, time_end, dtt

    ! Single parameter file from the command line (holds both domains).
    call yelmo_load_command_line_args(path_par)

    ! Shared timestepping (driver-owned) from the [ctrl] group.
    call nml_read(path_par, "ctrl", "tstep_method", tstep_method)
    call nml_read(path_par, "ctrl", "tstep_const",  tstep_const)
    call nml_read(path_par, "ctrl", "time_init",    time_init)
    call nml_read(path_par, "ctrl", "time_end",     time_end)
    call nml_read(path_par, "ctrl", "dtt",          dtt)
    call tstep_init(ts, time_init, time_end, method=tstep_method, units="year", &
                    time_ref=1950.0_wp, const_rel=tstep_const)

    ! Which hemispheres are active ([ctrl]).
    call nml_read(path_par, "ctrl", "active_north", active_north)
    call nml_read(path_par, "ctrl", "active_south", active_south)
    if (.not. (active_north .or. active_south)) then
        write(*,*) "yelmox_mgbi: error -- neither active_north nor active_south is set."
        stop 1
    end if

    ! === Shared, driver-owned barystatic sea level (one per run, common to both
    !     domains, exactly as in yelmox_bipolar). Restored once from a run-level
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
    if (active_north) call setup_domain(dom_north, "_north", outfldr_north)
    if (active_south) call setup_domain(dom_south, "_south", outfldr_south)

    ! === Output setup (shared cadence, from [tm_1D]/[tm_2D]) ===
    call timeout_init(tm_2D, path_par, "tm_2D", "heavy", time_init, time_end)
    call timeout_init(tm_1D, path_par, "tm_1D", "small", time_init, time_end)
    if (active_north) call write_domain_init(dom_north, outfldr_north)
    if (active_south) call write_domain_init(dom_south, outfldr_south)

    ! === Main time loop (shared timeline + dtt) ===
    call tstep_print_header(ts)
    do while (.not. ts%is_finished)
        call tstep_update(ts, dtt)
        call tstep_print(ts)

        ! Shared sea level: update once per step, before either domain advances,
        ! so both domains' isostasy see (and, under fastiso, contribute to) the
        ! same barystatic sea level.
        call bsl_update(bsl, ts%time_rel)

        if (active_north) call yelmox_step(dom_north, ts, bsl)
        if (active_south) call yelmox_step(dom_south, ts, bsl)

        ! Evaluate the shared output cadence once (timeout_check advances state).
        do_2D = tm_2D%active .and. timeout_check(tm_2D, ts%time)
        do_1D = tm_1D%active .and. timeout_check(tm_1D, ts%time)

        do_restart = .false.
        if (active_north) call write_domain_step(dom_north, outfldr_north)
        if (active_south) call write_domain_step(dom_south, outfldr_south)

        ! One shared bsl restart at the run root when either domain wrote a bundle.
        if (do_restart) &
            call bsl_restart_write(bsl, trim(restart_bundle_dir(ts%time))//"/bsl_restart.nc", ts%time)
    end do

    ! === Finalize: capture the final state + a final restart bundle per domain ===
    if (active_north) call write_domain_step(dom_north, outfldr_north, force=.true.)
    if (active_south) call write_domain_step(dom_south, outfldr_south, force=.true.)
    call bsl_restart_write(bsl, trim(restart_bundle_dir(ts%time))//"/bsl_restart.nc", ts%time)

    write(*,*)
    write(*,*) "yelmox_mgbi: run complete at time =", ts%time
    if (active_north) write(*,*) "  "//trim(dom_north%ctl%domain)//" H_ice max =", &
                                 maxval(dom_north%yelmo%tpo%now%H_ice)
    if (active_south) write(*,*) "  "//trim(dom_south%ctl%domain)//" H_ice max =", &
                                 maxval(dom_south%yelmo%tpo%now%H_ice)

    if (active_north) call yelmo_end(dom_north%yelmo, time=ts%time)
    if (active_south) call yelmo_end(dom_south%yelmo, time=ts%time)

contains

    subroutine setup_domain(dom, suffix, outfldr)
        ! Initialize one hemisphere: sub-models + hi-res hub + coupler maps, its
        ! output folder + regions, and the initial (cold) or restored state. The
        ! shared bsl was already initialized/restored by the driver above.
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: suffix
        character(len=*), intent(out)   :: outfldr

        call domain_init(dom, path_par, ts%time, group_suffix=suffix)

        ! Inject the driver-owned timeline values the domain logic needs.
        dom%ctl%tstep_method = tstep_method
        dom%ctl%dtt          = dtt

        ! Each domain writes to a subfolder named after its domain.
        outfldr = trim(dom%ctl%domain)//"/"
        call execute_command_line('mkdir -p "'//trim(outfldr)//'"')

        call domain_regions_init(dom, trim(outfldr))

        if (trim(dom%ctl%restart) == "None") then
            call domain_init_state(dom, ts, bsl)
        else
            call domain_restart_read(dom, trim(dom%ctl%restart), ts, bsl)
            call refresh_htopo(dom)
        end if

        write(*,*)
        write(*,*) "yelmox_mgbi: domain initialized ("//trim(adjustl(suffix))//")"
        write(*,*) "  domain      : "//trim(dom%ctl%domain)
        write(*,*) "  Yelmo grid  : "//trim(dom%ctl%grid_yelmo), dom%yelmo%grd%nx, dom%yelmo%grd%ny
        write(*,*) "  topo grid   : "//trim(dom%ctl%grid_name),  dom%topo%nx,      dom%topo%ny
        write(*,*) "  coupler maps: ", dom%cpl%nmaps
        write(*,*) "  output dir  : "//trim(outfldr)
        write(*,*)
    end subroutine setup_domain

    subroutine write_domain_init(dom, outfldr)
        ! Create output files and write the initial 2D + 1D records.
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr

        if (tm_2D%active) then
            call domain_write_init(dom, trim(outfldr), ts%time)
            call domain_write_step(dom, trim(outfldr), ts%time)
        end if
        if (tm_1D%active) call domain_write_1D(dom, trim(outfldr), ts%time, init=.TRUE.)
    end subroutine write_domain_init

    subroutine write_domain_step(dom, outfldr, force)
        ! Append 2D/1D records on the shared cadence (do_2D/do_1D), and write a
        ! restart bundle on the domain's dt_restart cadence -- flagging do_restart
        ! so the driver writes the single shared bsl restart. force=.true. writes
        ! everything unconditionally (final step).
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr
        logical, intent(in), optional   :: force

        logical :: do_force
        do_force = .false.
        if (present(force)) do_force = force

        if (do_2D .or. (do_force .and. tm_2D%active)) &
            call domain_write_step(dom, trim(outfldr), ts%time)
        if (do_1D .or. (do_force .and. tm_1D%active)) &
            call domain_write_1D(dom, trim(outfldr), ts%time)

        if (do_force) then
            call domain_restart_write(dom, ts%time, outfldr=trim(outfldr))
        else if (dom%ctl%dt_restart > 0.0_wp .and. &
                 mod(nint(ts%time*100), nint(dom%ctl%dt_restart*100)) == 0) then
            call domain_restart_write(dom, ts%time, outfldr=trim(outfldr))
            do_restart = .true.
        end if
    end subroutine write_domain_step

end program yelmox_mgbi
