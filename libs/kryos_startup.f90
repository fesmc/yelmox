module kryos_startup
    ! Initial state of a kryos_domain: the cold start (boundary state, Yelmo
    ! state init, cold-start ice state) and restart bundles (one folder per
    ! time holding one restart file per stateful component).

    use ncio
    use timestepping, only : tstep_class
    use yelmo,        only : wp, yelmo_update_equil, yelmo_init_state, yelmo_init_topo, &
                             yelmo_print_bound, yelmo_restart_write, yelmo_restart_init, &
                             yelmo_regions_update, mask_bed_ocean
    use yelmo_defs,   only : MASK_ICE_NONE
    use yelmo_tools,  only : smooth_gauss_2D
    use basal_dragging,   only : calc_cb_ref
    use yelmo_topography, only : calc_ytopo_diagnostic
    use yelmo_io,         only : yelmo_restart_read_topo_bnd
    use marine_shelf, only : marshelf_restart_write, marshelf_restart_read
    use fastisostasy, only : isos_init_ref, isos_init_state, isos_restart_write, &
                             bsl_class, bsl_update, bsl_restart_read, bsl_restart_write
    use smbpal,       only : smbpal_restart_write, smbpal_restart_read
    use kryos,        only : kryos_domain, remap, remap_method_smooth
    use kryos_regions,  only : domain_init_marine_ice
    use kryos_coupling, only : refresh_hub, step_climate, update_climate, step_smb, &
                               step_marine_shelf, couple_to_yelmo, couple_isostasy_to_yelmo, &
                               couple_smb_to_yelmo, check_isostasy_reference
    use kryos_forcing,  only : tsforcing_class, tsforcing_restart_write

    implicit none
    private

    public :: domain_startup, domain_init_state, domain_init_isostasy, domain_init_ice
    public :: domain_opt_init_cb_ref
    public :: bsl_startup, run_restart_write
    public :: domain_restart_write, domain_restart_read
    public :: restart_bundle_dir, restart_bundle_mkdir

contains

    subroutine bsl_startup(bsl, ts, fldr)
        ! Restore the shared, driver-owned barystatic sea level from a run-level
        ! restart bundle (fldr/bsl_restart.nc) and refresh it for the current
        ! time. No-op when fldr is "None" (bsl_init already set the cold state).
        ! bsl is prognostic under method="fastiso"/"mixed" and cannot be
        ! re-derived from time, hence the explicit restore.
        type(bsl_class),   intent(inout) :: bsl
        type(tstep_class), intent(in)    :: ts
        character(len=*),  intent(in)    :: fldr

        if (trim(fldr) == "None") return
        call bsl_restart_read(bsl, trim(fldr)//"/bsl_restart.nc")
        call bsl_update(bsl, ts%time_rel)
    end subroutine bsl_startup

    subroutine domain_startup(dom, ts, bsl, restore_bsl, tsf)
        ! Establish the domain state after domain_init: cold start (ctl%restart
        ! == "None") builds the initial boundary state; otherwise the restart
        ! bundle is restored and the hi-res hub rebuilt from the restored models.
        ! restore_bsl (default .true.) also restores the shared bsl from the same
        ! bundle folder -- the single-domain convention, where the run-level
        ! bsl_restart.nc lives in the domain's bundle. Multi-domain drivers
        ! restore the bsl once themselves (bsl_startup) and pass .false..
        ! Flavor drivers with their own cold start (esm, rembo) keep their own
        ! cold branch and call this for the restart branch only.
        ! The restart branch does not rebuild the climate/smb or marine-shelf
        ! forcing (not held in the bundle): every driver re-establishes it
        ! after this call with its own climate step + marine-shelf step.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        type(bsl_class),   intent(inout) :: bsl
        logical, intent(in), optional    :: restore_bsl
        type(tsforcing_class), intent(in), optional :: tsf   ! transient forcing (cold start)

        logical :: do_bsl

        do_bsl = .true.
        if (present(restore_bsl)) do_bsl = restore_bsl

        if (trim(dom%ctl%restart) == "None") then
            call domain_init_state(dom, ts, bsl, tsf=tsf)
        else
            if (do_bsl) call bsl_startup(bsl, ts, trim(dom%ctl%restart))
            call domain_restart_read(dom, trim(dom%ctl%restart), ts, bsl)
            call refresh_hub(dom)
        end if
    end subroutine domain_startup

    subroutine run_restart_write(dom, bsl, time, tsf, fldr)
        ! Single-domain restart: write the domain bundle + the run-level shared bsl
        ! restart (+ the driver-owned tsforcing state, when present) into one
        ! folder. By default the auto-named per-time folder; `fldr` overrides it
        ! (e.g. the forcing-increment "restart-<n>" folders). Multi-domain drivers
        ! write per-domain bundles + one run-root bsl bundle themselves.
        type(kryos_domain),      intent(inout)        :: dom
        type(bsl_class),       intent(inout)        :: bsl
        real(wp),              intent(in)           :: time
        type(tsforcing_class), intent(in), optional :: tsf
        character(len=*),      intent(in), optional :: fldr

        character(len=1024) :: bundle

        if (present(fldr)) then
            bundle = trim(fldr)
        else
            bundle = restart_bundle_dir(time)
        end if

        call domain_restart_write(dom, time, fldr=trim(bundle))
        call bsl_restart_write(bsl, trim(bundle)//"/bsl_restart.nc", time)
        if (present(tsf)) call tsforcing_restart_write(tsf, trim(bundle), time)
    end subroutine run_restart_write

    subroutine domain_init_state(dom, ts, bsl, tsf)
        ! Build the initial boundary state and initialize the Yelmo state
        ! variables, then run the domain-specific cold-start setup.
        ! bsl is the shared, driver-owned sea level; the driver has already called
        ! bsl_update for the initial time, so this routine only consumes it. The
        ! optional transient forcing (tsf) applies its initial anomalies to the
        ! startup climate, keeping the cold-start state consistent with the loop.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        type(bsl_class),   intent(inout) :: bsl
        type(tsforcing_class), intent(in), optional :: tsf

        ! Sea level + isostasy reference state (isostasy runs on grid_isos)
        call domain_init_isostasy(dom, ts, bsl)

        ! Refresh the hub from the initial geometry; climate/smb/mshlf read from it.
        call refresh_hub(dom)

        ! Climate on grid_clim (init: snapclim/snapesm at time_rel), then the
        ! surface mass balance on grid_smb (init=.true. runs the smbpal ITM
        ! equilibration before the first update).
        if (dom%ctl%with_climate) then
            call update_climate(dom, ts, tsf=tsf, init=.true.)
            call step_smb(dom, ts, init=.true.)
        end if

        ! Marine shelf through the (already refreshed) hub.
        call step_marine_shelf(dom, ts)

        ! Assemble the Yelmo boundary state from the freshly produced outputs.
        call couple_to_yelmo(dom)

        ! Cold-start friction guess for the optimization (restart restores cb_ref),
        ! set before the state init so its first dynamics solve already uses it.
        call domain_opt_init_cb_ref(dom)

        ! Initialize state variables (dyn, therm, mat) with a cold base
        call yelmo_print_bound(dom%yelmo%bnd)
        call yelmo_init_state(dom%yelmo, time=ts%time, thrm_method="robin-cold")

        ! Cold-start ice state ([coupling] init_marine_H, init_method). Cold
        ! start only; restart skips it.
        call domain_init_ice(dom, ts)

    end subroutine domain_init_state

    subroutine domain_opt_init_cb_ref(dom)
        ! Cold-start basal friction of the optimization (equil_method == "opt"; a
        ! restart restores cb_ref instead): opt.cf_init > 0 sets a uniform cb_ref,
        ! cf_init <= 0 starts from the till friction of the bed (cb_tgt, from the
        ! ytill parameters).
        type(kryos_domain), intent(inout) :: dom

        if (trim(dom%ctl%equil_method) /= "opt") return

        if (dom%opt%cf_init > 0.0_wp) then
            dom%yelmo%dyn%now%cb_ref = dom%opt%cf_init
        else
            call calc_cb_ref(dom%yelmo%dyn%now%cb_ref, dom%yelmo%bnd%z_bed, dom%yelmo%bnd%z_bed_sd, &
                    dom%yelmo%bnd%z_sl, dom%yelmo%bnd%H_sed, dom%yelmo%dyn%par%till_f_sed, &
                    dom%yelmo%dyn%par%till_sed_min, dom%yelmo%dyn%par%till_sed_max, &
                    dom%yelmo%dyn%par%till_cf_ref, dom%yelmo%dyn%par%till_cf_min, &
                    dom%yelmo%dyn%par%till_z0, dom%yelmo%dyn%par%till_z1, dom%yelmo%dyn%par%till_n_sd, &
                    dom%yelmo%dyn%par%till_scale_zb, dom%yelmo%dyn%par%till_scale_sed)
        end if
    end subroutine domain_opt_init_cb_ref

    subroutine domain_init_isostasy(dom, ts, bsl)
        ! Cold-start isostasy: reference and initial state on grid_isos from the
        ! Yelmo geometry, checked against Yelmo's own reference bedrock, then the
        ! bedrock / sea surface landed on the Yelmo grid. The ice load is
        ! coarsened conservatively (refined bilinearly), as in step_isostasy.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        type(bsl_class),   intent(inout) :: bsl

        real(wp), allocatable :: z_bed_ref_i(:,:), H_ice_ref_i(:,:)
        real(wp), allocatable :: z_bed_i(:,:), H_ice_i(:,:)
        character(len=256) :: gi, gy
        character(len=32)  :: mth_load

        gi = trim(dom%ctl%grid_isos)
        gy = trim(dom%ctl%grid_ice)

        mth_load = remap_method_smooth(real(dom%yelmo%grd%G%dx, wp), dom%ctl%dx_isos)
        call remap(dom, dom%yelmo%bnd%z_bed_ref, gy, z_bed_ref_i, gi, mth_load)
        call remap(dom, dom%yelmo%bnd%H_ice_ref, gy, H_ice_ref_i, gi, mth_load)
        call isos_init_ref(dom%isos, z_bed_ref_i, H_ice_ref_i)
        call remap(dom, dom%yelmo%bnd%z_bed,      gy, z_bed_i,     gi, mth_load)
        call remap(dom, dom%yelmo%tpo%now%H_ice,  gy, H_ice_i,     gi, mth_load)
        call isos_init_state(dom%isos, z_bed_i, H_ice_i, ts%time, bsl)
        call check_isostasy_reference(dom)
        call couple_isostasy_to_yelmo(dom)
    end subroutine domain_init_isostasy

    subroutine domain_init_ice(dom, ts)
        ! Cold-start ice state, after yelmo_init_state ([coupling]): no ice where
        ! the present-day bed is ocean (kill_shelves), LGM-like marine ice
        ! (init_marine_H), then init_method -- none, a short equilibration with
        ! constant boundaries (equil), or the ice reconstruction recon_path as the
        ! initial ice (recon) or only as the reference ice (recon_ref) -- and last
        ! an equilibration with topography fixed (time_equil_thrm > 0).
        type(kryos_domain), intent(inout) :: dom
        type(tstep_class),  intent(in)    :: ts

        if (dom%ctl%kill_shelves) then
            where (dom%yelmo%dta%pd%mask_bed == mask_bed_ocean) dom%yelmo%bnd%mask_ice = MASK_ICE_NONE
        end if

        if (dom%ctl%init_marine_H) call domain_init_marine_ice(dom)

        select case(trim(dom%ctl%init_method))
            case("none")
                ! Keep the initial state.
            case("equil")
                if (dom%ctl%with_ice_sheet) &
                    call yelmo_update_equil(dom%yelmo, ts%time, time_tot=dom%ctl%init_equil_time, &
                                            dt=1.0_wp, topo_fixed=.FALSE.)
            case("recon", "recon_ref")
                call domain_init_recon(dom, ts)
        end select

        if (dom%ctl%with_ice_sheet .and. dom%ctl%time_equil_thrm > 0.0_wp) &
            call yelmo_update_equil(dom%yelmo, ts%time, time_tot=dom%ctl%time_equil_thrm, &
                                    dt=dom%ctl%dtt, topo_fixed=.TRUE.)

    end subroutine domain_init_ice

    subroutine domain_init_recon(dom, ts)
        ! Start from an ice reconstruction (recon_path, e.g. ICE-6G_C at the LGM):
        ! it becomes the reference ice thickness and, with recon, the initial ice
        ! on the recon_codes regions (bed > -500 m). Then refresh the surface
        ! and (via the hub) the climate/smb, and stabilize the dynamic fields.
        type(kryos_domain), intent(inout) :: dom
        type(tstep_class),  intent(in)    :: ts

        integer  :: nx, ny, k
        logical, allocatable :: on_codes(:,:)
        real(wp) :: beta_min_save

        nx = dom%yelmo%tpo%par%nx
        ny = dom%yelmo%tpo%par%ny

        ! Load the reconstruction (slice 1) into the reference ice thickness.
        call nc_read(dom%ctl%recon_path, dom%ctl%recon_var, dom%yelmo%bnd%H_ice_ref, &
                     start=[1,1,1], count=[nx,ny,1])

        if (trim(dom%ctl%init_method) == "recon") then
            ! Initial ice from the reconstruction on the given regions.
            allocate(on_codes(nx,ny))
            on_codes = .false.
            do k = 1, dom%ctl%n_recon_codes
                where (abs(dom%yelmo%bnd%regions - dom%ctl%recon_codes(k)) < 1e-3_wp) on_codes = .true.
            end do
            where (dom%yelmo%bnd%z_bed > -500.0_wp .and. on_codes) &
                dom%yelmo%tpo%now%H_ice = dom%yelmo%bnd%H_ice_ref
            call smooth_gauss_2D(dom%yelmo%tpo%now%H_ice, dx=real(dom%yelmo%grd%G%dx,wp), f_sigma=2.0_wp)
            call yelmo_init_topo(dom%yelmo, trim(dom%ctl%path_par), &
                                 dom%yelmo%par%nml_init_topo, ts%time, load_topo=.FALSE.)
        end if

        ! Update surface topography fields (fixed H), then remove thin floating ice.
        call yelmo_update_equil(dom%yelmo, ts%time, time_tot=1.0_wp, dt=1.0_wp, topo_fixed=.TRUE.)
        where(dom%yelmo%tpo%now%mask_bed == 5 .and. dom%yelmo%tpo%now%H_ice < 50.0_wp) &
            dom%yelmo%tpo%now%H_ice = 0.0_wp
        call yelmo_update_equil(dom%yelmo, ts%time, time_tot=1.0_wp, dt=1.0_wp, topo_fixed=.TRUE.)

        if (trim(dom%ctl%init_method) == "recon") then
            ! Store the clean thickness as the reference state (drives smb masks).
            dom%yelmo%bnd%H_ice_ref = dom%yelmo%tpo%now%H_ice
        end if

        ! Refresh the hub and climate/smb to reflect the new geometry, then land
        ! the smb on the Yelmo grid for the stabilization below.
        call refresh_hub(dom)
        call step_climate(dom, ts)
        call step_smb(dom, ts)
        call couple_smb_to_yelmo(dom)

        ! Stabilize the dynamic fields with a raised beta_min.
        if (dom%ctl%with_ice_sheet) then
            beta_min_save = dom%yelmo%dyn%par%beta_min
            dom%yelmo%dyn%par%beta_min = 100.0_wp
            call yelmo_update_equil(dom%yelmo, ts%time, time_tot=2e2_wp, dt=5.0_wp, &
                                    topo_fixed=.FALSE.)
            dom%yelmo%dyn%par%beta_min = beta_min_save
        end if

    end subroutine domain_init_recon

    function restart_bundle_dir(time, outfldr) result(bundle)
        ! Auto-named per-time restart bundle folder: "<outfldr>restart-<kyr>-kyr".
        ! Shared by domain_restart_write and the driver (for the shared bsl bundle)
        ! so a domain's sub-model restarts and the run's bsl restart use identical
        ! folder naming.
        real(wp),         intent(in)           :: time
        character(len=*), intent(in), optional :: outfldr
        character(len=1024) :: bundle

        character(len=1024) :: prefix
        character(len=32)   :: time_str

        prefix = ""
        if (present(outfldr)) prefix = trim(outfldr)
        write(time_str,"(f20.3)") time*1e-3
        bundle = trim(prefix)//"restart-"//trim(adjustl(time_str))//"-kyr"
    end function restart_bundle_dir

    subroutine restart_bundle_mkdir(time, outfldr)
        ! Create the auto-named restart bundle folder (mkdir -p). The driver uses
        ! this for the shared bsl_restart.nc, which is written outside
        ! domain_restart_write (which creates its own per-domain bundle folder) and
        ! so needs its folder created explicitly.
        real(wp),         intent(in)           :: time
        character(len=*), intent(in), optional :: outfldr
        call execute_command_line('mkdir -p "'//trim(restart_bundle_dir(time, outfldr))//'"')
    end subroutine restart_bundle_mkdir

    subroutine domain_restart_write(dom, time, fldr, outfldr)
        ! Write a restart bundle: a folder (per time, or `fldr`) holding one
        ! restart file per stateful sub-model with fixed names. The hi-res hub is
        ! not written -- it is rebuilt by refresh_hub from the restored models.
        ! The shared barystatic sea level is NOT written here -- the driver owns it
        ! and writes a single bsl_restart.nc for the whole run.
        ! `outfldr` (optional) prefixes the auto-named per-time folder, so each
        ! domain of a multi-domain run writes into its own subfolder.
        type(kryos_domain), intent(inout) :: dom
        real(wp),         intent(in)    :: time
        character(len=*), intent(in), optional :: fldr
        character(len=*), intent(in), optional :: outfldr

        character(len=1024) :: bundle

        if (present(fldr)) then
            bundle = trim(fldr)
        else
            bundle = restart_bundle_dir(time, outfldr)
        end if

        call execute_command_line('mkdir -p "'//trim(bundle)//'"')

        call isos_restart_write(dom%isos,    trim(bundle)//"/isos_restart.nc",  time)
        ! Only checkpoint the ice-sheet state when the ice sheet is active. With
        ! with_ice_sheet=False the Yelmo dynamics never run, so there is no
        ! meaningful ice state to write (and the restart writer is not exercised
        ! in that mode).
        if (dom%ctl%with_ice_sheet) &
            call yelmo_restart_write(dom%yelmo,  trim(bundle)//"/yelmo_restart.nc", time)
        call marshelf_restart_write(dom%mshlf, trim(bundle)//"/marine_shelf.nc", time)
        call smbpal_restart_write(dom%smb,   trim(bundle)//"/smbpal_restart.nc", time)

        write(*,*) "domain_restart_write:: wrote bundle "//trim(bundle)
    end subroutine domain_restart_write

    subroutine domain_restart_read(dom, fldr, ts, bsl)
        ! Restore all stateful sub-models from a restart bundle folder. The shared
        ! barystatic sea level (bsl) is restored by the driver (bsl is prognostic
        ! under method="fastiso"/"mixed" and cannot be re-derived from time, so the
        ! driver reads bsl_now back from the run's bsl_restart.nc and calls
        ! bsl_update once); this routine only consumes the restored bsl.
        !
        ! Isostasy is restored through its proper init-from-restart path
        ! (isos_init_state with use_restart), NOT a bare isos_restart_read: the
        ! latter loads the state arrays but skips the post-read setup that
        ! isos_init_state performs (ODE state = now%w, calc_z_ss / calc_Haf /
        ! calc_masks, time_prognostics), without which the isostasy ODE solver
        ! restarts from an uninitialized state and the run is discontinuous.
        type(kryos_domain),  intent(inout) :: dom
        character(len=*),  intent(in)    :: fldr
        type(tstep_class), intent(in)    :: ts
        type(bsl_class),   intent(inout) :: bsl

        real(wp), allocatable :: z_bed_i(:,:), H_ice_i(:,:)
        character(len=256) :: gi, gy
        character(len=32)  :: mth_load

        gi = trim(dom%ctl%grid_isos)
        gy = trim(dom%ctl%grid_ice)

        ! Restore Yelmo first: it provides the current H_ice/z_bed for isostasy.
        ! Two steps, mirroring yelmo's native init-from-restart:
        !   - yelmo_restart_read_topo_bnd loads the geometry [tpo]+[bnd]
        !     (H_ice, z_bed, ...);
        !   - yelmo_restart_init loads [dyn,therm,mat] + mask_bed, activates the
        !     predictor-corrector and initializes the passive-tracer backends
        !     (elsa, tracer) from their sidecar files.
        ! The topo diagnostics (f_ice/f_grnd/H_grnd/z_srf) are reconciled below.
        !
        ! Only when the ice sheet is active. With with_ice_sheet=False the spin-up
        ! wrote no yelmo_restart.nc (see the matching guard in domain_restart_write),
        ! and the Yelmo dynamics never run; the geometry already loaded by
        ! domain_init (observed topography) is what the forcing-only run uses and
        ! what the isostasy/marine-shelf restore below consume.
        if (dom%ctl%with_ice_sheet) then
            call yelmo_restart_read_topo_bnd(dom%yelmo%tpo, dom%yelmo%bnd, dom%yelmo%time, &
                    dom%yelmo%par%restart_interpolated, dom%yelmo%grd, dom%yelmo%par%domain, &
                    dom%yelmo%par%grid_name, trim(fldr)//"/yelmo_restart.nc", ts%time)
            call yelmo_restart_init(dom%yelmo, trim(fldr)//"/yelmo_restart.nc", ts%time)
            dom%yelmo%par%use_restart = .true.
        end if

        ! Restore isostasy via isos_init_state (reads state + reference from the
        ! bundle and runs the full post-read setup), on the isos grid. The shared
        ! bsl was already restored + updated by the driver before this call.
        dom%isos%par%use_restart = .true.
        dom%isos%par%restart     = trim(fldr)//"/isos_restart.nc"
        mth_load = remap_method_smooth(real(dom%yelmo%grd%G%dx, wp), dom%ctl%dx_isos)
        call remap(dom, dom%yelmo%bnd%z_bed,     gy, z_bed_i, gi, mth_load)
        call remap(dom, dom%yelmo%tpo%now%H_ice, gy, H_ice_i, gi, mth_load)
        call isos_init_state(dom%isos, z_bed_i, H_ice_i, ts%time, bsl)
        call check_isostasy_reference(dom)
        call couple_isostasy_to_yelmo(dom)

        ! Restore marine shelf and the (prognostic, for ITM) snowpack state.
        call marshelf_restart_read(dom%mshlf, trim(fldr)//"/marine_shelf.nc")
        call smbpal_restart_read(dom%smb, trim(fldr)//"/smbpal_restart.nc")

        ! Reconcile Yelmo topo diagnostics (f_ice/f_grnd/H_grnd/z_srf) from the
        ! restored H_ice and the isostasy-updated z_bed/z_sl, then recompute the
        ! regional aggregates -- so the first 1D output after a restart reflects
        ! the restored state instead of the stale cold-start diagnostics.
        call calc_ytopo_diagnostic(dom%yelmo%tpo, dom%yelmo%dyn, dom%yelmo%mat, &
                                   dom%yelmo%thrm, dom%yelmo%bnd)
        call yelmo_regions_update(dom%yelmo)

        write(*,*) "domain_restart_read:: restored bundle "//trim(fldr)
    end subroutine domain_restart_read

end module kryos_startup
