program yelmox_bipolar
    ! Multigrid yelmox driver -- bipolar (two hemispheric domains + ocean box model).
    !
    ! Advances a Northern- and a Southern-Hemisphere ice_domain on a shared
    ! timeline, coupled through a shared Ocean Box Model (OBM), following the
    ! original yelmox_bipolar convention: one parameter file holds both domains,
    ! each domain's groups carry a hemisphere suffix (yelmo_north, coupling_south,
    ! snap_north, ...), while shared blocks ([ctrl], [barysealevel], the [obm]
    ! parameter block, the yelmo physics groups ydyn/ytopo/...) have no suffix.
    ! Distinct group names also let `runme -p group.name=val` target one
    ! hemisphere unambiguously.
    !
    ! [ctrl] active_north / active_south select which hemispheres run; a bipolar
    ! run is never more than these two domains, so they are held as two explicit
    ! ice_domain variables (not an array) -- the inter-domain ocean coupling is
    ! asymmetric (north <-> obm%fn/thetan/tn, south <-> obm%fs/thetas/ts).
    !
    ! Per-step coupling order (matches yelmox_bipolar): shared sea level, then per
    ! domain spinup + isostasy, then one OBM step, then per domain ice sheet +
    ! climate, then the ocean exchanges (atm->obm, ism->obm freshwater flux,
    ! hysteresis forcing, obm->ism ocean temperature), then per domain marine
    ! shelf (which reads the obm-updated ocean temperature). Per-domain physics
    ! live in libs/yelmox_domain.f90; the ocean coupling lives in
    ! yelmox_bipolar/obm_coupling.f90. See docs/multigrid.md.

    use nml
    use timestepping
    use timeout
    use yelmo,        only : yelmo_load_command_line_args, wp, yelmo_end
    use fastisostasy, only : bsl_class, bsl_init, bsl_update, bsl_restart_write
    use yelmox_domain
    use obm_defs,     only : obm_class
    use obm,          only : obm_init, obm_update, &
                             write_obm_init, write_obm_update, write_obm_restart
    use obm_coupling, only : obm_coupling_ctl, obm_ctl_load, obm_masks_init, obm_exchange

    implicit none

    character(len=512) :: path_par
    type(tstep_class)  :: ts
    type(bsl_class)    :: bsl                 ! shared, driver-owned sea level
    character(len=512) :: restart_bsl

    ! Two explicit hemispheric domains (bipolar: never more than north + south).
    type(ice_domain)   :: dom_north, dom_south
    character(len=512) :: outfldr_north, outfldr_south
    logical            :: active_north, active_south

    ! Shared, driver-owned ocean box model + its coupling control (obm_coupling).
    type(obm_class)        :: obm
    type(obm_coupling_ctl) :: oc
    character(len=512)     :: obm_file, obm_file_restart

    type(timeout_class) :: tm_2D, tm_2Dsm, tm_1D
    logical             :: do_2D, do_2Dsm, do_1D, do_restart
    logical             :: wrote_restart_north, wrote_restart_south

    real(wp) :: dtt

    ! Single parameter file from the command line (holds both domains + the obm).
    call yelmo_load_command_line_args(path_par)

    ! Shared timestepping (driver-owned) from the [ctrl] group.
    call timeline_init(ts, dtt, path_par, "ctrl")

    ! Which hemispheres are active ([ctrl]).
    call nml_read(path_par, "ctrl", "active_north", active_north)
    call nml_read(path_par, "ctrl", "active_south", active_south)
    if (.not. (active_north .or. active_south)) then
        write(*,*) "yelmox_bipolar: error -- neither active_north nor active_south is set."
        stop 1
    end if

    ! Ocean coupling control ([ctrl] switches + nautilus hysteresis parameters).
    call obm_ctl_load(oc, path_par)

    ! === Shared, driver-owned barystatic sea level (one per run, common to both
    !     domains, exactly as in yelmox_bipolar). Restored once from a run-level
    !     bsl_restart.nc ([ctrl] restart_bsl) when restarting. ===
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, ts%time_rel)
    call nml_read(path_par, "ctrl", "restart_bsl", restart_bsl)
    call bsl_startup(bsl, ts, restart_bsl)

    ! === Per-domain initialization ===
    if (active_north) call setup_domain(dom_north, "_north", outfldr_north)
    if (active_south) call setup_domain(dom_south, "_south", outfldr_south)

    ! Hydrographic masks (Yelmo grid) restricting the freshwater flux per domain.
    call obm_masks_init(oc, dom_north, dom_south, active_north, active_south)

    ! === Ocean box model init + first output record ===
    if (oc%active_obm) then
        obm_file         = trim(oc%obm_name)//".nc"
        obm_file_restart = trim(oc%obm_name)//"_restart.nc"
        call obm_init(obm, path_par, oc%obm_name)
        call write_obm_init(obm_file, ts%time, "years")
        call write_obm_update(obm, obm_file, oc%obm_name, ts%time)
    end if

    ! === Output setup (shared cadence, from [tm_2D]/[tm_2Dsm]/[tm_1D]) ===
    call timeout_init(tm_2D,   path_par, "tm_2D",   "heavy",  ts%time_init, ts%time_end)
    call timeout_init(tm_2Dsm, path_par, "tm_2Dsm", "medium", ts%time_init, ts%time_end)
    call timeout_init(tm_1D,   path_par, "tm_1D",   "small",  ts%time_init, ts%time_end)
    if (active_north) call write_domain_init(dom_north, outfldr_north)
    if (active_south) call write_domain_init(dom_south, outfldr_south)

    ! === Main time loop (shared timeline + dtt) ===
    call tstep_print_header(ts)
    do while (.not. ts%is_finished)
        call tstep_update(ts, dtt)
        call tstep_print(ts)

        ! Shared sea level: update once per step, before either domain advances.
        call bsl_update(bsl, ts%time_rel)

        ! Spinup relaxation + isostasy (both domains) -- before the OBM step.
        if (active_north) call advance_isostasy(dom_north)
        if (active_south) call advance_isostasy(dom_south)

        ! Ocean box model: one step, using last step's freshwater/atmos forcing.
        if (oc%active_obm) call obm_update(obm, dtt, oc%obm_name)

        ! Ice sheet + hi-res hub refresh + climate/smb (both domains).
        if (active_north) call advance_dynamics(dom_north)
        if (active_south) call advance_dynamics(dom_south)

        ! Inter-domain ocean coupling (shared obm): atm->obm, ism->obm freshwater
        ! flux, hysteresis forcing, obm->ism ocean temperature.
        call obm_exchange(oc, obm, dom_north, dom_south, active_north, active_south, &
                          ts%time, ts%time_init, dtt)

        ! Marine shelf (both domains) -- reads the obm-updated snapclim to_ann.
        if (active_north) call step_marine_shelf(dom_north, ts)
        if (active_south) call step_marine_shelf(dom_south, ts)

        ! === Output (shared cadence; timeout_check advances state, call once) ===
        do_2D   = tm_2D%active   .and. timeout_check(tm_2D, ts%time)
        do_2Dsm = tm_2Dsm%active .and. timeout_check(tm_2Dsm, ts%time)
        do_1D   = tm_1D%active   .and. timeout_check(tm_1D, ts%time)

        wrote_restart_north = .false.
        wrote_restart_south = .false.
        if (active_north) call write_domain_step(dom_north, outfldr_north, wrote_restart_north)
        if (active_south) call write_domain_step(dom_south, outfldr_south, wrote_restart_south)
        if (oc%active_obm .and. do_1D) call write_obm_update(obm, obm_file, oc%obm_name, ts%time)

        ! Shared bsl (+ obm) restart at the run root when either domain wrote one.
        do_restart = wrote_restart_north .or. wrote_restart_south
        if (do_restart) then
            call restart_bundle_mkdir(ts%time)
            call bsl_restart_write(bsl, trim(restart_bundle_dir(ts%time))//"/bsl_restart.nc", ts%time)
            if (oc%active_obm) call write_obm_restart(obm, obm_file_restart, ts%time, "years")
        end if
    end do

    ! === Finalize: capture the final state + a final restart bundle per domain ===
    if (active_north) call write_domain_step(dom_north, outfldr_north, wrote_restart_north, force=.true.)
    if (active_south) call write_domain_step(dom_south, outfldr_south, wrote_restart_south, force=.true.)
    call restart_bundle_mkdir(ts%time)
    call bsl_restart_write(bsl, trim(restart_bundle_dir(ts%time))//"/bsl_restart.nc", ts%time)
    if (oc%active_obm) call write_obm_restart(obm, obm_file_restart, ts%time, "years")

    write(*,*)
    write(*,*) "yelmox_bipolar: run complete at time =", ts%time
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

        ! Each domain writes to a subfolder named after its domain.
        outfldr = trim(dom%ctl%domain)//"/"
        call execute_command_line('mkdir -p "'//trim(outfldr)//'"')

        call domain_regions_init(dom, trim(outfldr))

        ! Cold start or per-domain restart; the shared bsl was already
        ! initialized/restored once by the driver (bsl_startup above).
        call domain_startup(dom, ts, bsl, restore_bsl=.false.)

        write(*,*)
        write(*,*) "yelmox_bipolar: domain initialized ("//trim(adjustl(suffix))//")"
        write(*,*) "  domain      : "//trim(dom%ctl%domain)
        write(*,*) "  Yelmo grid  : "//trim(dom%ctl%grid_yelmo), dom%yelmo%grd%G%nx, dom%yelmo%grd%G%ny
        write(*,*) "  topo grid   : "//trim(dom%ctl%grid_name),  dom%topo%nx,      dom%topo%ny
        write(*,*) "  coupler maps: ", dom%cpl%nmaps
        write(*,*) "  output dir  : "//trim(outfldr)
        write(*,*)
    end subroutine setup_domain

    subroutine advance_isostasy(dom)
        ! Per-domain part that precedes the OBM step: spinup relaxation +
        ! cb_ref/tf_corr tuning, then isostasy against the shared sea level.
        type(ice_domain), intent(inout) :: dom
        call step_optimize(dom, ts)
        call step_isostasy(dom, ts, bsl)
    end subroutine advance_isostasy

    subroutine advance_dynamics(dom)
        ! Per-domain part after the OBM step and before the ocean coupling: ice
        ! sheet update, hi-res hub refresh, then climate + surface mass balance.
        type(ice_domain), intent(inout) :: dom
        call step_icesheet(dom, ts)
        call refresh_htopo(dom)
        call step_climate(dom, ts)
    end subroutine advance_dynamics

    subroutine write_domain_init(dom, outfldr)
        ! Create output files and write the initial 2D + 1D records.
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr

        if (tm_2D%active) then
            call domain_write_init(dom, trim(outfldr), ts%time)
            call domain_write_step(dom, trim(outfldr), ts%time)
        end if
        if (tm_2Dsm%active) then
            call domain_write_init_sm(dom, trim(outfldr), ts%time)
            call domain_write_step_sm(dom, trim(outfldr), ts%time)
        end if
        if (tm_1D%active) call domain_write_1D(dom, trim(outfldr), ts%time, init=.TRUE.)
    end subroutine write_domain_init

    subroutine write_domain_step(dom, outfldr, wrote_restart, force)
        ! Append 2D/1D records on the shared cadence (do_2D/do_1D), and write a
        ! restart bundle on the domain's dt_restart cadence -- reporting it via
        ! wrote_restart so the driver writes the single shared bsl (+ obm)
        ! restart. force=.true. writes everything unconditionally (final step).
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr
        logical,          intent(out)   :: wrote_restart
        logical, intent(in), optional   :: force

        logical :: do_force
        do_force = .false.
        if (present(force)) do_force = force

        if (do_2D .or. (do_force .and. tm_2D%active)) &
            call domain_write_step(dom, trim(outfldr), ts%time)
        if (do_2Dsm .or. (do_force .and. tm_2Dsm%active)) &
            call domain_write_step_sm(dom, trim(outfldr), ts%time)
        if (do_1D .or. (do_force .and. tm_1D%active)) &
            call domain_write_1D(dom, trim(outfldr), ts%time)

        wrote_restart = .false.
        if (do_force) then
            call domain_restart_write(dom, ts%time, outfldr=trim(outfldr))
        else if (tstep_due(ts%time, dom%ctl%dt_restart)) then
            call domain_restart_write(dom, ts%time, outfldr=trim(outfldr))
            wrote_restart = .true.
        end if
    end subroutine write_domain_step

end program yelmox_bipolar
