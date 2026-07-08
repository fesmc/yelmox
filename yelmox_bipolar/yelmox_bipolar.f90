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
    use ncio,         only : nc_read
    use timestepping
    use timeout
    use yelmo,        only : yelmo_load_command_line_args, wp, yelmo_end
    use fastisostasy, only : bsl_class, bsl_init, bsl_update, &
                             bsl_restart_read, bsl_restart_write
    use yelmox_domain
    use obm_defs,     only : obm_class
    use obm,          only : obm_init, obm_update, &
                             write_obm_init, write_obm_update, write_obm_restart
    use obm_coupling, only : coupling_atm2obm, coupling_ism2obm, coupling_obm2ism, &
                             update_bipolar_hyster_forcing

    implicit none

    character(len=512) :: path_par
    type(tstep_class)  :: ts
    type(bsl_class)    :: bsl                 ! shared, driver-owned sea level
    character(len=512) :: restart_bsl

    ! Two explicit hemispheric domains (bipolar: never more than north + south).
    type(ice_domain)   :: dom_north, dom_south
    character(len=512) :: outfldr_north, outfldr_south
    logical            :: active_north, active_south

    ! Shared, driver-owned ocean box model + its coupling switches.
    type(obm_class)    :: obm
    logical            :: active_obm, ism2obm, obm2ism, atm2obm
    character(len=512) :: obm_name, obm_file, obm_file_restart
    logical            :: couple_fwf_north, couple_fwf_south
    character(len=512) :: fwf_definition
    character(len=512) :: hydro_mask_north_path, hydro_mask_south_path
    real(wp), allocatable :: hydro_mask_north(:,:), hydro_mask_south(:,:)
    logical            :: hyster_on
    character(len=512) :: hyster_forcing, hyster_forcing_method
    real(wp)           :: hyster_rate, hyster_positive_branch_time

    type(timeout_class) :: tm_2D, tm_1D
    logical             :: do_2D, do_1D, do_restart

    character(len=56) :: tstep_method
    real(wp)          :: tstep_const, time_init, time_end, dtt

    ! Single parameter file from the command line (holds both domains + the obm).
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
        write(*,*) "yelmox_bipolar: error -- neither active_north nor active_south is set."
        stop 1
    end if

    ! Ocean coupling switches ([ctrl]).
    call read_obm_control()

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

    ! Hydrographic masks (Yelmo grid) restricting the freshwater flux per domain.
    if (ism2obm) then
        if (couple_fwf_north .and. active_north) then
            allocate(hydro_mask_north(dom_north%yelmo%grd%G%nx, dom_north%yelmo%grd%G%ny))
            call nc_read(hydro_mask_north_path, "mask", hydro_mask_north)
        end if
        if (couple_fwf_south .and. active_south) then
            allocate(hydro_mask_south(dom_south%yelmo%grd%G%nx, dom_south%yelmo%grd%G%ny))
            call nc_read(hydro_mask_south_path, "mask", hydro_mask_south)
        end if
    end if

    ! === Ocean box model init + first output record ===
    if (active_obm) then
        obm_file         = trim(obm_name)//".nc"
        obm_file_restart = trim(obm_name)//"_restart.nc"
        call obm_init(obm, path_par, obm_name)
        call write_obm_init(obm_file, ts%time, "years")
        call write_obm_update(obm, obm_file, obm_name, ts%time)
    end if

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

        ! Shared sea level: update once per step, before either domain advances.
        call bsl_update(bsl, ts%time_rel)

        ! Spinup relaxation + isostasy (both domains) -- before the OBM step.
        if (active_north) call advance_isostasy(dom_north)
        if (active_south) call advance_isostasy(dom_south)

        ! Ocean box model: one step, using last step's freshwater/atmos forcing.
        if (active_obm) call obm_update(obm, dtt, obm_name)

        ! Ice sheet + hi-res hub refresh + climate/smb (both domains).
        if (active_north) call advance_dynamics(dom_north)
        if (active_south) call advance_dynamics(dom_south)

        ! === Inter-domain ocean coupling (shared obm) ===
        if (atm2obm) then
            if (active_north) call coupling_atm2obm(dom_north, obm, "north", ts%time)
            if (active_south) call coupling_atm2obm(dom_south, obm, "south", ts%time)
        end if
        if (ism2obm) then
            if (couple_fwf_north .and. active_north) &
                call coupling_ism2obm(dom_north, obm, hydro_mask_north, "north", fwf_definition)
            if (couple_fwf_south .and. active_south) &
                call coupling_ism2obm(dom_south, obm, hydro_mask_south, "south", fwf_definition)
        end if
        if (trim(obm_name) == "nautilus" .and. hyster_on) &
            call update_bipolar_hyster_forcing(ts%time, time_init, obm, dtt, &
                    hyster_positive_branch_time, hyster_rate, hyster_forcing, hyster_forcing_method)
        if (obm2ism) then
            if (active_north) call coupling_obm2ism(dom_north, obm, obm_name, "north")
            if (active_south) call coupling_obm2ism(dom_south, obm, obm_name, "south")
        end if

        ! Marine shelf (both domains) -- reads the obm-updated snapclim to_ann.
        if (active_north) call step_marine_shelf(dom_north, ts)
        if (active_south) call step_marine_shelf(dom_south, ts)

        ! === Output (shared cadence; timeout_check advances state, call once) ===
        do_2D = tm_2D%active .and. timeout_check(tm_2D, ts%time)
        do_1D = tm_1D%active .and. timeout_check(tm_1D, ts%time)

        do_restart = .false.
        if (active_north) call write_domain_step(dom_north, outfldr_north)
        if (active_south) call write_domain_step(dom_south, outfldr_south)
        if (active_obm .and. do_1D) call write_obm_update(obm, obm_file, obm_name, ts%time)

        ! Shared bsl (+ obm) restart at the run root when either domain wrote one.
        if (do_restart) then
            call restart_bundle_mkdir(ts%time)
            call bsl_restart_write(bsl, trim(restart_bundle_dir(ts%time))//"/bsl_restart.nc", ts%time)
            if (active_obm) call write_obm_restart(obm, obm_file_restart, ts%time, "years")
        end if
    end do

    ! === Finalize: capture the final state + a final restart bundle per domain ===
    if (active_north) call write_domain_step(dom_north, outfldr_north, force=.true.)
    if (active_south) call write_domain_step(dom_south, outfldr_south, force=.true.)
    call restart_bundle_mkdir(ts%time)
    call bsl_restart_write(bsl, trim(restart_bundle_dir(ts%time))//"/bsl_restart.nc", ts%time)
    if (active_obm) call write_obm_restart(obm, obm_file_restart, ts%time, "years")

    write(*,*)
    write(*,*) "yelmox_bipolar: run complete at time =", ts%time
    if (active_north) write(*,*) "  "//trim(dom_north%ctl%domain)//" H_ice max =", &
                                 maxval(dom_north%yelmo%tpo%now%H_ice)
    if (active_south) write(*,*) "  "//trim(dom_south%ctl%domain)//" H_ice max =", &
                                 maxval(dom_south%yelmo%tpo%now%H_ice)

    if (active_north) call yelmo_end(dom_north%yelmo, time=ts%time)
    if (active_south) call yelmo_end(dom_south%yelmo, time=ts%time)

contains

    subroutine read_obm_control()
        ! Load the ocean-coupling switches ([ctrl]) and, when needed, the
        ! freshwater-flux masks and nautilus hysteresis-forcing parameters.
        active_obm = .false.
        ism2obm    = .false.
        obm2ism    = .false.
        atm2obm    = .false.
        obm_name   = "none"
        call nml_read(path_par, "ctrl", "active_obm", active_obm)
        call nml_read(path_par, "ctrl", "ism2obm",    ism2obm)
        call nml_read(path_par, "ctrl", "obm2ism",    obm2ism)
        call nml_read(path_par, "ctrl", "atm2obm",    atm2obm)
        call nml_read(path_par, "ctrl", "obm_name",   obm_name)

        couple_fwf_north = .false.
        couple_fwf_south = .false.
        fwf_definition   = "dVdt"
        if (ism2obm) then
            call nml_read(path_par, "ctrl", "couple_fwf_north", couple_fwf_north)
            call nml_read(path_par, "ctrl", "couple_fwf_south", couple_fwf_south)
            call nml_read(path_par, "ctrl", "fwf_definition",   fwf_definition)
            if (couple_fwf_north) &
                call nml_read(path_par, "ctrl", "hydro_mask_north", hydro_mask_north_path)
            if (couple_fwf_south) &
                call nml_read(path_par, "ctrl", "hydro_mask_south", hydro_mask_south_path)
        end if

        hyster_on = .false.
        if (trim(obm_name) == "nautilus") then
            call nml_read(path_par, obm_name, "hyster_on",                   hyster_on)
            call nml_read(path_par, obm_name, "hyster_forcing",              hyster_forcing)
            call nml_read(path_par, obm_name, "hyster_forcing_method",       hyster_forcing_method)
            call nml_read(path_par, obm_name, "hyster_rate",                 hyster_rate)
            call nml_read(path_par, obm_name, "hyster_positive_branch_time", hyster_positive_branch_time)
        end if
    end subroutine read_obm_control

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
        if (tm_1D%active) call domain_write_1D(dom, trim(outfldr), ts%time, init=.TRUE.)
    end subroutine write_domain_init

    subroutine write_domain_step(dom, outfldr, force)
        ! Append 2D/1D records on the shared cadence (do_2D/do_1D), and write a
        ! restart bundle on the domain's dt_restart cadence -- flagging do_restart
        ! so the driver writes the single shared bsl (+ obm) restart. force=.true.
        ! writes everything unconditionally (final step).
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

end program yelmox_bipolar
