program yelmox_bipolar
    ! Multigrid yelmox driver -- bipolar (two hemispheric domains + ocean box model).
    !
    ! Advances a Northern- and a Southern-Hemisphere kryos_domain on a shared
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
    ! kryos_domain variables (not an array) -- the inter-domain ocean coupling is
    ! asymmetric (north <-> obm%fn/thetan/tn, south <-> obm%fs/thetas/ts).
    !
    ! Per-step coupling order: shared sea level, then per domain spin-up tuning +
    ! isostasy, then one OBM step, then per domain ice sheet + climate + smb,
    ! then the ocean exchanges (atm->obm, ism->obm freshwater flux, hysteresis
    ! forcing, obm->ism ocean temperature), then per domain marine shelf (which
    ! reads the obm-updated ocean temperature). Per-domain physics live in the
    ! kryos modules (libs/kryos*.f90); the ocean coupling lives in
    ! yelmox_bipolar/obm_coupling.f90. See docs/multigrid.md.

    use nml
    use timestepping
    use timeout
    use yelmo,        only : yelmo_load_command_line_args, wp, yelmo_end
    use fastisostasy, only : bsl_class, bsl_init, bsl_update, bsl_restart_write
    use kryos,          only : kryos_domain, domain_init
    use kryos_regions,  only : domain_regions_init
    use kryos_coupling, only : step_spinup_tuning, step_isostasy, couple_to_yelmo, &
                               step_icesheet, couple_yelmo_to_htopo, step_climate, &
                               step_smb, step_marine_shelf
    use kryos_startup,  only : domain_startup, bsl_startup, domain_restart_write, &
                               restart_bundle_dir, restart_bundle_mkdir
    use kryos_output,   only : domain_write_init, domain_write_step, &
                               domain_write_init_sm, domain_write_step_sm, &
                               domain_write_1D
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
    type(kryos_domain)   :: dom_north, dom_south
    character(len=512) :: outfldr_north, outfldr_south
    logical            :: active_north, active_south

    ! Shared, driver-owned ocean box model + its coupling control (obm_coupling).
    type(obm_class)        :: obox      ! ocean box model state ("obm" is the module name)
    type(obm_coupling_ctl) :: oc
    character(len=512)     :: obm_file, obm_restart

    type(timeout_class) :: tm_2D, tm_2Dsm, tm_1D, tm_rst
    logical             :: do_2D, do_2Dsm, do_1D, do_rst

    real(wp) :: dtt

    ! Single parameter file from the command line (holds both domains + the obm).
    call yelmo_load_command_line_args(path_par)

    ! Shared timestepping (driver-owned) from the [ctrl] group.
    call tstep_init(ts, path_par, "ctrl", dtt)

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
    !     domains). Restored once from the run-level restart bundle ([ctrl]
    !     restart_bsl, which also holds the OBM restart) when restarting. ===
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, ts%time_rel)
    call nml_read(path_par, "ctrl", "restart_bsl", restart_bsl)
    call bsl_startup(bsl, ts, restart_bsl)

    ! === Per-domain initialization ===
    if (active_north) call setup_domain(dom_north, "_north", outfldr_north)
    if (active_south) call setup_domain(dom_south, "_south", outfldr_south)

    ! Hydrographic masks (Yelmo grid) restricting the freshwater flux per domain.
    call obm_masks_init(oc, dom_north, dom_south, active_north, active_south)

    ! === Ocean box model init (from the run-level bundle when restarting) +
    !     output file ===
    if (oc%active_obm) then
        obm_file    = trim(oc%obm_name)//".nc"
        obm_restart = "None"
        if (trim(restart_bsl) /= "None") obm_restart = trim(restart_bsl)//"/obm_restart.nc"
        call obm_init(obox, path_par, oc%obm_name, obm_restart)
        call write_obm_init(obm_file, ts%time, "years")
    end if

    ! === Output + restart schedules (shared cadence, from [tm_2D]/[tm_2Dsm]/
    !     [tm_1D]/[tm_rst]; a restart bundle is always written at time_end) ===
    call timeout_init(tm_2D,   path_par, "tm_2D",   "heavy",   ts%time_init, ts%time_end)
    call timeout_init(tm_2Dsm, path_par, "tm_2Dsm", "medium",  ts%time_init, ts%time_end)
    call timeout_init(tm_1D,   path_par, "tm_1D",   "small",   ts%time_init, ts%time_end)
    call timeout_init(tm_rst,  path_par, "tm_rst",  "restart", ts%time_init, ts%time_end)
    if (active_north) call write_domain_init(dom_north, outfldr_north)
    if (active_south) call write_domain_init(dom_south, outfldr_south)

    ! === Main time loop (shared timeline + dtt) ===
    ! Output and restarts are written at the top of the loop for the current
    ! time (time_init on the first pass), then the loop exits once the run is
    ! finished, else the time is advanced and both domains stepped. The final
    ! state and a final restart bundle are always written.
    call tstep_print_header(ts)
    do

        ! === Output + restarts (shared cadence) ===
        do_2D   = tm_2D%active   .and. (timeout_check(tm_2D,   ts%time) .or. ts%is_finished)
        do_2Dsm = tm_2Dsm%active .and. (timeout_check(tm_2Dsm, ts%time) .or. ts%is_finished)
        do_1D   = tm_1D%active   .and. (timeout_check(tm_1D,   ts%time) .or. ts%is_finished)
        do_rst  = timeout_check(tm_rst, ts%time) .or. ts%is_finished

        if (active_north) call write_domain_step(dom_north, outfldr_north)
        if (active_south) call write_domain_step(dom_south, outfldr_south)
        if (oc%active_obm .and. do_1D) call write_obm_update(obox, obm_file, oc%obm_name, ts%time)

        ! Run-level restart bundle at the run root (shared bsl + obm), next to
        ! the per-domain bundles.
        if (do_rst) then
            call restart_bundle_mkdir(ts%time)
            call bsl_restart_write(bsl, trim(restart_bundle_dir(ts%time))//"/bsl_restart.nc", ts%time)
            if (oc%active_obm) call write_obm_restart(obox, &
                    trim(restart_bundle_dir(ts%time))//"/obm_restart.nc", ts%time, "years")
        end if

        if (ts%is_finished) exit

        call tstep_update(ts, dtt)
        call tstep_print(ts)

        ! Shared sea level: update once per step, before either domain advances.
        call bsl_update(bsl, ts%time_rel)

        ! Spin-up tuning + isostasy (each domain), before the OBM step.
        if (active_north) then
            call step_spinup_tuning(dom_north, ts)
            call step_isostasy(dom_north, ts, bsl)
        end if
        if (active_south) then
            call step_spinup_tuning(dom_south, ts)
            call step_isostasy(dom_south, ts, bsl)
        end if

        ! Ocean box model: one step, using last step's freshwater/atmos forcing.
        if (oc%active_obm) call obm_update(obox, dtt, oc%obm_name)

        ! Ice sheet, hub refresh, climate + surface mass balance (each domain).
        if (active_north) then
            call couple_to_yelmo(dom_north)   ! bedrock now; smb + shelf melt lag one step
            call step_icesheet(dom_north, ts)
            call couple_yelmo_to_htopo(dom_north)
            call step_climate(dom_north, ts)
            call step_smb(dom_north, ts)
        end if
        if (active_south) then
            call couple_to_yelmo(dom_south)
            call step_icesheet(dom_south, ts)
            call couple_yelmo_to_htopo(dom_south)
            call step_climate(dom_south, ts)
            call step_smb(dom_south, ts)
        end if

        ! Inter-domain ocean coupling (shared obm): atm->obm, ism->obm freshwater
        ! flux, hysteresis forcing, obm->ism ocean temperature.
        call obm_exchange(oc, obox, dom_north, dom_south, active_north, active_south, &
                          ts%time, ts%time_init, dtt)

        ! Marine shelf (both domains) -- reads the obm-updated snapclim to_ann.
        if (active_north) call step_marine_shelf(dom_north, ts)
        if (active_south) call step_marine_shelf(dom_south, ts)
    end do

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
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: suffix
        character(len=*), intent(out)   :: outfldr

        call domain_init(dom, path_par, ts%time, group_suffix=suffix)

        ! Each domain writes to a subfolder named after its domain.
        outfldr = trim(dom%ctl%domain)//"/"
        call execute_command_line('mkdir -p "'//trim(outfldr)//'"')

        call domain_regions_init(dom, trim(outfldr))

        ! Cold start or per-domain restart; the shared bsl was already
        ! initialized/restored once by the driver (bsl_startup above). On a
        ! restart, re-establish the climate/smb and marine-shelf forcing from the
        ! restored state (the bundle does not hold them), so the first step and
        ! the first output see a valid boundary state.
        call domain_startup(dom, ts, bsl, restore_bsl=.false.)
        if (trim(dom%ctl%restart) /= "None") then
            call step_climate(dom, ts)
            call step_smb(dom, ts)
            call step_marine_shelf(dom, ts)
        end if

        write(*,*)
        write(*,*) "yelmox_bipolar: domain initialized ("//trim(adjustl(suffix))//")"
        write(*,*) "  domain      : "//trim(dom%ctl%domain)
        write(*,*) "  Yelmo grid  : "//trim(dom%ctl%grid_ice), dom%yelmo%grd%G%nx, dom%yelmo%grd%G%ny
        write(*,*) "  topo grid   : "//trim(dom%ctl%grid_hub),  dom%topo%nx,      dom%topo%ny
        write(*,*) "  coupler maps: ", dom%cpl%nmaps
        write(*,*) "  output dir  : "//trim(outfldr)
        write(*,*)
    end subroutine setup_domain

    subroutine write_domain_init(dom, outfldr)
        ! Create the 2D + 1D output files.
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr

        if (tm_2D%active)   call domain_write_init(dom, trim(outfldr), ts%time)
        if (tm_2Dsm%active) call domain_write_init_sm(dom, trim(outfldr), ts%time)
        if (tm_1D%active)   call domain_write_1D(dom, trim(outfldr), ts%time, init=.TRUE.)
    end subroutine write_domain_init

    subroutine write_domain_step(dom, outfldr)
        ! Append 2D/1D records and write the domain restart bundle on the shared
        ! cadence (do_2D/do_2Dsm/do_1D/do_rst, set by the driver at the top of
        ! the loop); the driver writes the single shared bsl (+ obm) restart.
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr

        if (do_2D)   call domain_write_step(dom, trim(outfldr), ts%time)
        if (do_2Dsm) call domain_write_step_sm(dom, trim(outfldr), ts%time)
        if (do_1D)   call domain_write_1D(dom, trim(outfldr), ts%time)
        if (do_rst)  call domain_restart_write(dom, ts%time, outfldr=trim(outfldr))
    end subroutine write_domain_step

end program yelmox_bipolar
