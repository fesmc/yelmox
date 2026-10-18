program yelmox_bipolar
    ! Multigrid yelmox driver -- bipolar (two hemispheric domains + ocean box model).
    !
    ! Advances a Northern- and a Southern-Hemisphere kryos_domain on a shared
    ! timeline, coupled through a shared Ocean Box Model (OBM), following the
    ! original yelmox_bipolar convention: one parameter file holds both domains,
    ! each domain's groups carry a hemisphere suffix (yelmo_north, comps_south,
    ! snap_north, ...), while shared blocks ([ctrl], [barysealevel], the [obm]
    ! parameter block, the yelmo physics groups ydyn/ytopo/...) have no suffix.
    ! Distinct group names also let `runme -p group.name=val` target one
    ! hemisphere unambiguously.
    !
    ! [ctrl] active_north / active_south select which hemispheres run (none: the
    ! ocean box model alone, active_obm); a bipolar
    ! run is never more than these two domains, so they are held as two explicit
    ! kryos_domain variables (not an array) -- the inter-domain ocean coupling is
    ! asymmetric (north <-> obm%fn/thetan/tn, south <-> obm%fs/thetas/ts).
    !
    ! Per-step coupling order: shared sea level and indices ([indices], read
    ! by the climate of both domains and the OBM forcing), then per domain
    ! relaxation, optimization + isostasy, then one OBM step, then per domain
    ! ice sheet + climate + smb, then the ocean exchanges (atm->obm from the
    ! global index of [obm_forcing], ism->obm freshwater flux, hysteresis
    ! forcing, obm->ism ocean temperature through each domain's ocn_obm<sfx>),
    ! then per domain marine shelf (which reads the obm-updated ocean
    ! temperature). Per-domain physics live in the
    ! kryos modules (libs/kryos*.f90); the ocean coupling lives in
    ! yelmox_bipolar/obm_coupling.f90. See docs/multigrid.md.

    use nml
    use timestepping
    use timeout
    use precision,    only : dp
    use yelmo,        only : yelmo_load_command_line_args, wp, yelmo_end
    use fastisostasy, only : bsl_class, bsl_init, bsl_update, bsl_restart_write
    use kryos,          only : kryos_domain, domain_init
    use yelmox_climate, only : climate_end
    use kryos_regions,  only : domain_regions_init
    use kryos_coupling, only : step_relax, step_optimize, step_isostasy, couple_to_yelmo, &
                               step_icesheet, couple_yelmo_to_htopo, step_climate, &
                               step_surface, step_shelf
    use kryos_startup,  only : domain_startup, bsl_startup, domain_restart_write, &
                               restart_bundle_dir
    use kryos_output,   only : domain_write_init, domain_write_step, &
                               domain_write_init_sm, domain_write_step_sm, &
                               domain_write_1D
    use kryos_forcing,  only : indices_class, indices_init, indices_update, indices_kill, &
                               indices_restart_due, indices_restart_fldr, &
                               indices_restart_read, indices_restart_write, indices_write_step
    use obm_defs,     only : obm_class
    use obm,          only : obm_init, obm_update, &
                             write_obm_init, write_obm_update, write_obm_restart
    use obm_coupling, only : obm_coupling_ctl, obm_ctl_load, obm_coupling_init, obm_ocean_init, &
                             obm_exchange
    use timer,        only : timer_class, timer_init, timer_step, &
                             timer_print_totals, timer_write_totals

    implicit none

    character(len=512) :: path_par
    type(tstep_class)  :: ts
    type(bsl_class)    :: bsl                 ! shared, driver-owned sea level
    type(indices_class) :: idx                ! shared, driver-owned indices ([indices])
    real(wp)           :: fvar                ! response variable of the indices: ice volume [Gt]
    character(len=512) :: restart             ! run restart bundle ([ctrl]), or "None"

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

    ! Wall time of each part of the time loop (fesm-utils timer, both domains
    ! together), accumulated over the run: printed at the end, and written to
    ! timing.txt (tm_1D cadence).
    type(timer_class) :: tmr
    integer, parameter :: tt_output = 1, tt_forcing = 2, tt_relax = 3, tt_optimize = 4, &
                          tt_isostasy = 5, tt_obm = 6, tt_couple_yelmo = 7, tt_icesheet = 8, &
                          tt_couple_htopo = 9, tt_climate = 10, tt_surface = 11, &
                          tt_obm_exchange = 12, tt_shelf = 13
    logical :: timing_init = .TRUE.   ! timing.txt is created on the first 1D write

    real(dp) :: dtt

    ! Single parameter file from the command line (holds both domains + the obm).
    call yelmo_load_command_line_args(path_par)

    ! Shared timestepping (driver-owned) from the [ctrl] group.
    call tstep_init(ts, path_par, "ctrl", dtt)

    ! Which hemispheres are active ([ctrl]).
    call nml_read(path_par, "ctrl", "active_north", active_north)
    call nml_read(path_par, "ctrl", "active_south", active_south)

    ! Ocean coupling control ([ctrl] switches, [obm_forcing], nautilus
    ! hysteresis parameters); checked by obm_coupling_init below.
    call obm_ctl_load(oc, path_par)

    ! At least one component runs: an ice-sheet domain or the ocean box model
    ! alone.
    if (.not. (active_north .or. active_south .or. oc%active_obm)) then
        write(*,*) "yelmox_bipolar: error -- none of active_north, active_south, active_obm is set."
        stop 1
    end if

    ! === Shared, driver-owned barystatic sea level (one per run, common to both
    !     domains). Restored once from the run restart bundle ([ctrl] restart,
    !     which also holds the OBM restart and a subfolder per domain). ===
    call bsl_init(bsl, path_par, ts%time_rel)
    call bsl_update(bsl, ts%time_rel)
    call nml_read(path_par, "ctrl", "restart", restart)
    call bsl_startup(bsl, restart)

    ! === Shared, driver-owned indices on the forcing clock (both domains read
    !     them), restored from the run restart bundle. Before the domains: their
    !     climate checks the indices it names. ===
    call indices_init(idx, path_par, ts%time_rel)
    if (trim(restart) /= "None") call indices_restart_read(idx, trim(restart))

    ! === Per-domain initialization ===
    if (active_north) call setup_domain(dom_north, "_north", outfldr_north)
    if (active_south) call setup_domain(dom_south, "_south", outfldr_south)

    ! The ocean couplings checked against the indices and the domains
    ! (ocn_obm<sfx>), and the hydrographic masks (Yelmo grid) restricting the
    ! freshwater flux per domain.
    call obm_coupling_init(oc, idx, dom_north, dom_south, active_north, active_south)

    ! === Ocean box model init (from the run-level bundle when restarting) +
    !     output file ===
    if (oc%active_obm) then
        obm_file    = trim(oc%obm_name)//".nc"
        obm_restart = "None"
        if (trim(restart) /= "None") obm_restart = trim(restart)//"/obm_restart.nc"
        call obm_init(obox, path_par, oc%obm_name, obm_restart)
        call write_obm_init(obm_file, ts%time, "years")

        ! The domains' ocean from the box (ocn_obm<sfx>) and their shelf
        ! forcing, so the first record and the first step see it.
        call obm_ocean_init(oc, obox, dom_north, dom_south, active_north, active_south, ts)
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
    call timer_init(tmr, [character(len=16) :: "output", "forcing", "relax", "optimize", &
                          "isostasy", "obm", "couple_yelmo", "icesheet", "couple_htopo", &
                          "climate", "surface", "obm_exchange", "shelf"])

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
        if (do_1D) then
            call timer_write_totals(tmr, ts%time, "s", "timing.txt", init=timing_init)
            timing_init = .FALSE.
        end if

        ! Run restart bundle (auto-named per-time folder).
        if (do_rst) call write_run_restart(restart_bundle_dir(ts%time))

        call timer_step(tmr, tt_output)

        if (ts%is_finished) exit

        call tstep_update(ts, dtt)
        call tstep_print(ts)

        ! Shared sea level and indices: update once per step, before either
        ! domain advances. Response variable of the indices: the ice volume of
        ! the active domains [Gt].
        call bsl_update(bsl, ts%time_rel)
        fvar = 0.0_wp
        if (active_north) fvar = fvar + dom_north%yelmo%reg%V_ice * dom_north%yelmo%bnd%c%rho_ice * 1e-3_wp
        if (active_south) fvar = fvar + dom_south%yelmo%reg%V_ice * dom_south%yelmo%bnd%c%rho_ice * 1e-3_wp
        call indices_update(idx, ts%time_rel, var=fvar)
        call timer_step(tmr, tt_forcing)

        ! Relaxation + optimization + isostasy (each domain), before the OBM step.
        if (active_north) then
            call step_relax(dom_north, ts)
            call timer_step(tmr, tt_relax)
            call step_optimize(dom_north, ts)
            call timer_step(tmr, tt_optimize)
            call step_isostasy(dom_north, ts, bsl)
            call timer_step(tmr, tt_isostasy)
        end if
        if (active_south) then
            call step_relax(dom_south, ts)
            call timer_step(tmr, tt_relax)
            call step_optimize(dom_south, ts)
            call timer_step(tmr, tt_optimize)
            call step_isostasy(dom_south, ts, bsl)
            call timer_step(tmr, tt_isostasy)
        end if

        ! Ocean box model: one step, using last step's freshwater/atmos forcing.
        if (oc%active_obm) call obm_update(obox, real(dtt,wp), oc%obm_name)
        call timer_step(tmr, tt_obm)

        ! Ice sheet, hub refresh, climate + surface mass balance (each domain).
        if (active_north) then
            call couple_to_yelmo(dom_north)   ! bedrock now; smb + shelf melt lag one step
            call timer_step(tmr, tt_couple_yelmo)
            call step_icesheet(dom_north, ts)
            call timer_step(tmr, tt_icesheet)
            call couple_yelmo_to_htopo(dom_north)
            call timer_step(tmr, tt_couple_htopo)
            call step_climate(dom_north, ts, idx)
            call timer_step(tmr, tt_climate)
            call step_surface(dom_north, ts)
            call timer_step(tmr, tt_surface)
        end if
        if (active_south) then
            call couple_to_yelmo(dom_south)
            call timer_step(tmr, tt_couple_yelmo)
            call step_icesheet(dom_south, ts)
            call timer_step(tmr, tt_icesheet)
            call couple_yelmo_to_htopo(dom_south)
            call timer_step(tmr, tt_couple_htopo)
            call step_climate(dom_south, ts, idx)
            call timer_step(tmr, tt_climate)
            call step_surface(dom_south, ts)
            call timer_step(tmr, tt_surface)
        end if

        ! Inter-domain ocean coupling (shared obm): atm->obm (global index),
        ! ism->obm freshwater flux, hysteresis forcing, obm->ism ocean
        ! temperature (ocn_obm<sfx>).
        call obm_exchange(oc, obox, idx, dom_north, dom_south, active_north, active_south, &
                          ts%time, ts%time_init, dtt)
        call timer_step(tmr, tt_obm_exchange)

        ! Marine shelf (both domains) -- reads the obm-updated to_ann (ocn_obm<sfx>).
        if (active_north) call step_shelf(dom_north, ts)
        if (active_south) call step_shelf(dom_south, ts)
        call timer_step(tmr, tt_shelf)

        ! Forcing-increment restart each time an index moves by more than its
        ! restart_every_df (folders restart-<n>).
        if (indices_restart_due(idx)) then
            write(*,*) "yelmox_bipolar: forcing-increment restart "// &
                       trim(indices_restart_fldr(idx))//" at time =", ts%time
            call write_run_restart(indices_restart_fldr(idx))
        end if

        ! Kill switch: stop once the response has equilibrated at a forcing
        ! bound (any index). The top of the loop then writes the final state and exits.
        if (indices_kill(idx)) then
            write(*,*) "yelmox_bipolar: index kill switch tripped at time =", ts%time
            ts%is_finished = .TRUE.
        end if
    end do

    write(*,*)
    write(*,*) "yelmox_bipolar: run complete at time =", ts%time
    if (active_north) write(*,*) "  "//trim(dom_north%ctl%domain)//" H_ice max =", &
                                 maxval(dom_north%yelmo%tpo%now%H_ice)
    if (active_south) write(*,*) "  "//trim(dom_south%ctl%domain)//" H_ice max =", &
                                 maxval(dom_south%yelmo%tpo%now%H_ice)
    write(*,*)
    write(*,*) "yelmox_bipolar: wall time of the time loop by component"
    call timer_print_totals(tmr, "s")

    if (active_north) call yelmo_end(dom_north%yelmo, time=ts%time)
    if (active_south) call yelmo_end(dom_south%yelmo, time=ts%time)
    if (active_north) call climate_end(dom_north%cl)
    if (active_south) call climate_end(dom_south%cl)

contains

    subroutine setup_domain(dom, suffix, outfldr)
        ! Initialize one hemisphere: sub-models + hi-res hub + coupler maps, its
        ! output folder + regions, and the initial (cold) or restored state. The
        ! shared bsl was already initialized/restored by the driver above.
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: suffix
        character(len=*), intent(out)   :: outfldr

        call domain_init(dom, path_par, ts%time, idx, group_suffix=suffix)

        ! Each domain writes to a subfolder named after its domain.
        outfldr = trim(dom%ctl%domain)//"/"
        call execute_command_line('mkdir -p "'//trim(outfldr)//'"')

        call domain_regions_init(dom, trim(outfldr))

        ! Cold start or restart from the domain's subfolder of the run bundle;
        ! the shared bsl was already initialized/restored once by the driver
        ! (bsl_startup above). On a restart, re-establish the climate/smb and
        ! marine-shelf forcing from the restored state (the bundle does not hold
        ! them), so the first step and the first output see a valid boundary state.
        if (trim(restart) == "None") then
            call domain_startup(dom, ts, bsl, "None", idx, restore_bsl=.false.)
        else
            call domain_startup(dom, ts, bsl, trim(restart)//"/"//trim(dom%ctl%domain), idx, &
                                restore_bsl=.false.)
            call step_climate(dom, ts, idx)
            call step_surface(dom, ts)
            call step_shelf(dom, ts)
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
        ! Append 2D/1D records (the shared indices in each domain's 1D file) on
        ! the shared cadence (do_2D/do_2Dsm/do_1D, set by the driver at the top
        ! of the loop).
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr

        if (do_2D)   call domain_write_step(dom, trim(outfldr), ts%time)
        if (do_2Dsm) call domain_write_step_sm(dom, trim(outfldr), ts%time)
        if (do_1D)   call domain_write_1D(dom, trim(outfldr), ts%time)
        if (do_1D)   call indices_write_step(idx, dom%yelmo%reg%fnm, ts%time)
    end subroutine write_domain_step

    subroutine write_run_restart(bundle)
        ! Run restart bundle in the folder `bundle`: the shared bsl, indices and
        ! obm, and each active domain in its subfolder named after the domain.
        character(len=*), intent(in) :: bundle

        call execute_command_line('mkdir -p "'//trim(bundle)//'"')
        call bsl_restart_write(bsl, trim(bundle)//"/bsl_restart.nc", ts%time)
        call indices_restart_write(idx, trim(bundle), ts%time)
        if (oc%active_obm) call write_obm_restart(obox, trim(bundle)//"/obm_restart.nc", &
                                                  ts%time, "years")
        if (active_north) call domain_restart_write(dom_north, ts%time, &
                                fldr=trim(bundle)//"/"//trim(dom_north%ctl%domain))
        if (active_south) call domain_restart_write(dom_south, ts%time, &
                                fldr=trim(bundle)//"/"//trim(dom_south%ctl%domain))
    end subroutine write_run_restart

end program yelmox_bipolar
