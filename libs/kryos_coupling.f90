module kryos_coupling
    ! The per-step coupling primitives of a kryos_domain. Each step_* advances
    ! one component on its own grid, reading geometry from the hi-res hub; each
    ! couple_*_to_yelmo lands one component's output on the Yelmo grid, before
    ! Yelmo runs. Every primitive is a no-op when its component is inactive.

    use timestepping, only : tstep_class
    use precision,    only : dp
    use yelmo,        only : wp, yelmo_update
    use marine_shelf, only : marshelf_update, marshelf_update_shelf
    use fastisostasy, only : isos_update, bsl_class
    use yelmox_climate, only : climate_update
    use smbpal,       only : smbpal_update_monthly, smbpal_update_monthly_equil
    use surface_chion, only : surface_chion_spinup, surface_chion_update
    use smb_simple_m, only : smb_simple_set_mask, smb_simple_update
    use htopo,        only : htopo_update
    use ice_optimization, only : relax_update, optimize_cb_ref, optimize_tf_corr, &
                                 optimize_tf_corr_basin
    use kryos,        only : kryos_domain, remap, remap_method_smooth, cadence_due
    use kryos_regions, only : negis_update_cb_ref, calc_glacial_smb
    use kryos_forcing, only : coefs_class
    use ice_bnd,      only : ice_bnd_update

    implicit none
    private

    public :: step_relax, step_optimize, step_isostasy, step_ice_bnd, step_icesheet, step_climate, step_shelf
    public :: step_surface, couple_yelmo_to_htopo
    public :: couple_to_yelmo
    public :: couple_isostasy_to_yelmo, couple_surface_to_yelmo, couple_shelf_to_yelmo
    public :: couple_climate_to_yelmo
    public :: check_isostasy_reference

contains

    subroutine step_relax(dom, ts)
        ! Relaxation of the topography ([sim] relax): ytopo.topo_rel is
        ! relax.topo_rel with the timescale ramp tau1 -> tau2 until relax.time2
        ! (time elapsed since the start of the run), then the par-file values.
        type(kryos_domain), intent(inout) :: dom
        type(tstep_class),  intent(in)    :: ts

        if (.not. dom%ctl%relax) return

        call relax_update(dom%rlx%par, dom%yelmo%tpo%par, ts%time_elapsed, &
                          dom%rlx%topo_rel0, dom%rlx%topo_rel_tau0)
    end subroutine step_relax

    subroutine step_optimize(dom, ts, dtt)
        ! Optimization ([sim] opt), within the [opt] time windows: nudge the
        ! basal-friction field cb_ref and the marine thermal-forcing correction
        ! tf_corr toward present-day observations, over the step dtt.
        !
        ! cb_ref is a Yelmo-grid control, optimized in place. tf_corr is a
        ! marine_shelf control, optimized in place on grid_shelf, where it lives
        ! and is never remapped (a round trip through the Yelmo grid would smooth
        ! it every step). The Yelmo-grid inputs of its optimization (the current
        ! and observed geometry) are remapped onto grid_shelf instead; the
        ! basins are marine_shelf's own. At identity grids these remaps are copies.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        real(dp),          intent(in)    :: dtt   ! [yr] the driver's time step

        real(wp), allocatable :: H_ice_m(:,:), H_grnd_m(:,:), dHidt_m(:,:)
        real(wp), allocatable :: H_ice_pd_m(:,:), H_grnd_pd_m(:,:)
        character(len=256) :: gm, gy
        character(len=32)  :: mth

        if (.not. dom%ctl%opt) return

        gm = trim(dom%ctl%grid_shelf)
        gy = trim(dom%ctl%grid_ice)

        ! Basal friction (cb_ref) optimization -- Yelmo grid, in place.
        if (trim(dom%opt%opt_cf) == "L21" .and. ts%time_elapsed >= dom%opt%cf_time_init &
                            .and. ts%time_elapsed <= dom%opt%cf_time_end) then
            call optimize_cb_ref(dom%yelmo%dyn%now%cb_ref, dom%yelmo%tpo%now%H_ice, &
                    dom%yelmo%tpo%now%dHidt, dom%yelmo%bnd%z_bed, dom%yelmo%bnd%z_sl, &
                    dom%yelmo%dyn%now%ux_s, dom%yelmo%dyn%now%uy_s, &
                    dom%yelmo%dta%pd%H_ice, dom%yelmo%dta%pd%uxy_s, dom%yelmo%dta%pd%H_grnd, &
                    dom%opt%cf_min, dom%opt%cf_max, dom%yelmo%tpo%par%dx, &
                    dom%opt%sigma_err, dom%opt%sigma_vel, dom%opt%tau_c, dom%opt%H0, &
                    dt=dtt, fill_method=dom%opt%fill_method, fill_dist=dom%opt%sigma_err, &
                    cb_tgt=dom%yelmo%dyn%now%cb_tgt)
        end if

        ! Thermal-forcing correction (tf_corr) optimization -- grid_shelf, in
        ! place on mshlf (incl. its restart). The Yelmo-grid inputs land on
        ! grid_shelf averaged over the target cell when the Yelmo grid is finer
        ! (conservative), else interpolated (remap_method_smooth). The basins are
        ! marine_shelf's (mshlf%now%basins, on grid_shelf): the same regions layer
        ! and ids as Yelmo's bnd%basins when [marine_shelf] basins = [yelmo_masks]
        ! basins, as in every par file, so the ids of tf_basins.
        ! Method: "L21" (one correction per basin) or "L21-points".
        if (trim(dom%opt%opt_tf) /= "none" .and. ts%time_elapsed >= dom%opt%tf_time_init &
                            .and. ts%time_elapsed <= dom%opt%tf_time_end) then
            mth = remap_method_smooth(real(dom%yelmo%grd%G%dx, wp), dom%ctl%dx_shelf)
            call remap(dom, dom%yelmo%tpo%now%H_ice,  gy, H_ice_m,    gm, mth)
            call remap(dom, dom%yelmo%tpo%now%H_grnd, gy, H_grnd_m,   gm, mth)
            call remap(dom, dom%yelmo%tpo%now%dHidt,  gy, dHidt_m,    gm, mth)
            call remap(dom, dom%yelmo%dta%pd%H_ice,   gy, H_ice_pd_m, gm, mth)
            select case (trim(dom%opt%opt_tf))
            case ("L21")
                call optimize_tf_corr_basin(dom%mshlf%now%tf_corr, H_ice_m, H_grnd_m, &
                        dHidt_m, H_ice_pd_m, dom%mshlf%now%basins, &
                        dom%opt%H_grnd_lim, dom%opt%tau_m, dom%opt%m_temp, dom%opt%tf_min, dom%opt%tf_max, &
                        dom%opt%tf_basins, dt=dtt)
            case ("L21-points")
                call remap(dom, dom%yelmo%dta%pd%H_grnd, gy, H_grnd_pd_m, gm, mth)
                call optimize_tf_corr(dom%mshlf%now%tf_corr, H_ice_m, H_grnd_m, &
                        dHidt_m, H_ice_pd_m, H_grnd_pd_m, &
                        dom%opt%H_grnd_lim, dom%mshlf%now%basins, dom%opt%basin_fill, &
                        dom%opt%tau_m, dom%opt%m_temp, dom%opt%tf_min, dom%opt%tf_max, &
                        dom%ctl%dx_shelf, sigma=dom%opt%tf_sigma, dt=dtt)
            end select
        end if
    end subroutine step_optimize

    subroutine step_isostasy(dom, ts, bsl)
        ! Run isostasy on its own grid, on the dt_isos cadence (FastIsostasy
        ! integrates the whole interval; the bedrock is held in between): ice load
        ! from Yelmo, conservative onto a coarser grid_isos, bilinear onto a finer
        ! one (remap_method_smooth; a copy when equal). The bedrock /
        ! sea-surface outputs stay on grid_isos (in dom%isos%out); they are landed
        ! on the Yelmo grid by couple_isostasy_to_yelmo (in couple_to_yelmo, before
        ! step_icesheet). bsl is the shared, driver-owned sea level (already
        ! updated for this step by the driver); isos_update reads it and, under
        ! fastiso/mixed, writes back the prognostic bsl_now -- so with several
        ! domains sharing one bsl the sea level integrates every domain's ice load.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        type(bsl_class),   intent(inout) :: bsl

        real(wp), allocatable :: H_ice_i(:,:), dwdt_i(:,:)
        character(len=256) :: gi, gy
        character(len=32)  :: mth_load

        if (.not. dom%ctl%with_isostasy) return
        if (.not. cadence_due(ts%time_elapsed, dom%ctl%dt_isos)) return

        gi = trim(dom%ctl%grid_isos)
        gy = trim(dom%ctl%grid_ice)

        ! ice load + correction: Yelmo -> isos grid
        mth_load = remap_method_smooth(real(dom%yelmo%grd%G%dx, wp), dom%ctl%dx_isos)
        call remap(dom, dom%yelmo%tpo%now%H_ice,  gy, H_ice_i, gi, mth_load)
        call remap(dom, dom%yelmo%bnd%dzbdt_corr, gy, dwdt_i,  gi, mth_load)

        call isos_update(dom%isos, H_ice_i, ts%time, bsl, dwdt_corr=dwdt_i)
    end subroutine step_isostasy

    ! --- Yelmo-input couplers -------------------------------------------------
    ! Each coupler remaps one module's output onto the Yelmo grid and assigns it
    ! into yelmo%bnd, i.e. "remap what Yelmo needs, before Yelmo runs". They are
    ! called together by couple_to_yelmo (in the per-step sequence, before
    ! step_icesheet) and from the init/restart paths (before yelmo_init_state),
    ! so the Yelmo boundary assembly lives in one place. Each is a no-op when its
    ! component is inactive. At identity grids the remaps are copies.

    subroutine couple_to_yelmo(dom)
        ! Assemble the Yelmo boundary state from every coupled component. In the
        ! time loop isostasy was produced this step; smb and marine shelf were
        ! produced last step (the one-step coupling lag).
        type(kryos_domain), intent(inout) :: dom

        call couple_isostasy_to_yelmo(dom)
        call couple_surface_to_yelmo(dom)
        call couple_shelf_to_yelmo(dom)
        call couple_climate_to_yelmo(dom)
    end subroutine couple_to_yelmo

    subroutine couple_climate_to_yelmo(dom)
        ! Yelmo inputs the climate supplies directly: subglacial discharge (Qd,
        ! e.g. Greenland frontal melt, grid_clim -> Yelmo, conservative), when the
        ! backend provides it.
        type(kryos_domain), intent(inout) :: dom

        real(wp), allocatable :: Qd_y(:,:)

        if (.not. dom%clim%has_Qd) return

        call remap(dom, dom%clim%now%Qd, dom%ctl%grid_clim, Qd_y, dom%ctl%grid_ice, "con")
        dom%yelmo%bnd%Qd = Qd_y
    end subroutine couple_climate_to_yelmo

    subroutine couple_isostasy_to_yelmo(dom)
        ! Bedrock + sea surface from isostasy (grid_isos -> Yelmo). Only the
        ! isostatic *anomaly* crosses grids: Yelmo keeps the reference bedrock it
        ! read natively on its own grid and adds the remapped displacement,
        !
        !     z_bed = z_bed_ref + (w + we)
        !     z_sl  = bsl       + (z_ss - bsl)
        !
        ! rather than taking isostasy's absolute z_bed, which is built internally as
        ! ref%z_bed + w + we on a reference that was itself coarsened onto grid_isos.
        ! Passing the absolute field would replace Yelmo's bedrock with a
        ! grid_isos-resolution copy of it, smoothing away every trough, sill and
        ! pinning point finer than grid_isos -- at t = 0 (w = we = 0) that is a pure
        ! round trip carrying no isostatic signal at all. The displacement and the
        ! sea-surface perturbation are long-wavelength, and are the only part of the
        ! solution that is genuinely insensitive to the grid it was solved on.
        !
        ! Yelmo applies the same decomposition when restarting from an interpolated
        ! file (yelmo_ice.f90, "isostatic offset from z_bed_ref"), and FastIsostasy
        ! likewise refuses to round-trip z_bed through the ice grid on restart.
        !
        ! For z_sl the reference is the scalar bsl, which is uniform and so needs no
        ! remapping; the spatial part z_ss - bsl (= isos ref%z_ss + dz_ss) is what
        ! crosses grids.
        type(kryos_domain), intent(inout) :: dom

        real(wp), allocatable :: dz_bed_y(:,:), dz_ss_y(:,:)
        character(len=256) :: gi, gy
        character(len=32)  :: mth

        if (.not. dom%ctl%with_isostasy) return

        gi  = trim(dom%ctl%grid_isos)
        gy  = trim(dom%ctl%grid_ice)
        mth = remap_method_smooth(dom%ctl%dx_isos, real(dom%yelmo%grd%G%dx, wp))

        call remap(dom, dom%isos%out%w + dom%isos%out%we,   gi, dz_bed_y, gy, mth)
        call remap(dom, dom%isos%out%z_ss - dom%isos%now%bsl, gi, dz_ss_y, gy, mth)

        dom%yelmo%bnd%z_bed = dom%yelmo%bnd%z_bed_ref + dz_bed_y
        dom%yelmo%bnd%z_sl  = dom%isos%now%bsl        + dz_ss_y
    end subroutine couple_isostasy_to_yelmo

    subroutine check_isostasy_reference(dom)
        ! Verify that the reference bedrock isostasy is working from is the same
        ! field Yelmo holds natively, because couple_isostasy_to_yelmo rebuilds
        ! z_bed = z_bed_ref + (w + we) from the two of them independently.
        !
        ! They are paired by construction on a cold start (isos_init_ref is handed a
        ! remapped yelmo%bnd%z_bed_ref), but on a restart isostasy restores its own
        ! reference from isos_restart.nc, which may have come from a different
        ! topography. Nothing else would catch that: the run would proceed with a
        ! displacement field measured against one bedrock and applied to another.
        !
        ! The comparison is made on the isostasy grid: Yelmo's z_bed_ref is remapped
        ! there exactly as the reference was built (domain_init_isostasy), and
        ! isostasy's reference is out%z_bed - (w + we). The same bedrock then agrees
        ! to round-off on any isostasy grid; a different one differs by metres or more.
        type(kryos_domain), intent(inout) :: dom

        real(wp), parameter :: tol = 1.0_wp   ! [m] max |difference|: round-off only

        real(wp), allocatable :: z_bed_ref_i(:,:)
        character(len=256) :: gi, gy
        character(len=32)  :: mth
        real(wp) :: dmean, dmax

        if (.not. dom%ctl%with_isostasy) return

        gi  = trim(dom%ctl%grid_isos)
        gy  = trim(dom%ctl%grid_ice)
        mth = remap_method_smooth(real(dom%yelmo%grd%G%dx, wp), dom%ctl%dx_isos)

        call remap(dom, dom%yelmo%bnd%z_bed_ref, gy, z_bed_ref_i, gi, mth)

        dmean = sum(dom%isos%out%z_bed - dom%isos%out%w - dom%isos%out%we - z_bed_ref_i) &
                / real(size(z_bed_ref_i), wp)
        dmax  = maxval(abs(dom%isos%out%z_bed - dom%isos%out%w - dom%isos%out%we - z_bed_ref_i))

        write(*,*) "check_isostasy_reference:: z_bed_ref (isos - yelmo), on the isostasy grid [m]"
        write(*,*) "    grids:      ", trim(gy), " -> ", trim(gi), " (", trim(mth), ")"
        write(*,*) "    mean diff:  ", dmean
        write(*,*) "    max |diff|: ", dmax
        write(*,*) "    tolerance:  ", tol

        if (dmax > tol) then
            write(*,*) ""
            write(*,*) "check_isostasy_reference:: error: the isostasy reference bedrock does &
                       &not match yelmo%bnd%z_bed_ref."
            write(*,*) "  couple_isostasy_to_yelmo adds the isostatic displacement to Yelmo's own &
                       &reference, so the two must describe the same bedrock."
            write(*,*) "  On a restart this usually means isos_restart.nc came from a run with a &
                       &different topography than the one Yelmo is reading now."
            error stop 1
        end if
    end subroutine check_isostasy_reference

    subroutine couple_surface_to_yelmo(dom)
        ! Surface mass balance + surface temperature from the active SMB model
        ! (grid_surface -> Yelmo, conservative), with the we->ie unit scaling and the
        ! optional Greenland modifications. The producing step (step_surface,
        ! leaves smb/tsrf on grid_surface in the SMB model's own fields, or the
        ! climate has them (surface_method = "climate"); this coupler is the
        ! single place that lands them on Yelmo.
        type(kryos_domain), intent(inout) :: dom

        real(wp), allocatable :: smb_y(:,:), tsrf_y(:,:), ta_y(:,:), ta_pd_y(:,:)
        real(wp), allocatable :: smb_s(:,:), tsrf_s(:,:)
        character(len=256) :: gs, gc, gy

        if (.not. dom%ctl%with_surface) return

        gs = trim(dom%ctl%grid_surface)
        gc = trim(dom%ctl%grid_clim)
        gy = trim(dom%ctl%grid_ice)

        select case(trim(dom%ctl%surface_method))
            case("climate")
                ! The climate's own surface mass balance and surface temperature,
                ! through grid_surface as the surface models' fields.
                call remap(dom, dom%clim%now%smb,  gc, smb_s,  gs, "bilin")
                call remap(dom, dom%clim%now%tsrf, gc, tsrf_s, gs, "bilin")
                call remap(dom, smb_s,  gs, smb_y,  gy, "con")
                call remap(dom, tsrf_s, gs, tsrf_y, gy, "con")
            case("smb_simple")
                call remap(dom, dom%smbs%smb,   gs, smb_y,  gy, "con")
                call remap(dom, dom%smbs%t_srf, gs, tsrf_y, gy, "con")
            case("chion")
                call remap(dom, dom%schn%ann%smb,  gs, smb_y,  gy, "con")
                call remap(dom, dom%schn%ann%tsrf, gs, tsrf_y, gy, "con")
            case("smbpal")
                call remap(dom, dom%smb%ann%smb,  gs, smb_y,  gy, "con")
                call remap(dom, dom%smb%ann%tsrf, gs, tsrf_y, gy, "con")
        end select

        dom%yelmo%bnd%smb   = smb_y * dom%yelmo%bnd%c%conv_we_ie * 1e-3
        dom%yelmo%bnd%T_srf = tsrf_y

        ! Glacial-smb modification: reduce large negative smb toward a
        ! quasi glacial-interglacial index. Operates on the aggregated Yelmo-grid smb.
        if (dom%ctl%scale_glacial_smb) then
            call remap(dom, dom%clim%now%ta_ann, gc, ta_y,    gy, "bilin")
            call remap(dom, dom%clim%ref%ta_ann, gc, ta_pd_y, gy, "bilin")
            call calc_glacial_smb(dom%yelmo%bnd%smb, real(dom%yelmo%grd%lat,wp), ta_y, ta_pd_y, dom%gsmb)
        end if

        ! Limit to present-day ice extent: impose extra melt (4 m ie/a) wherever
        ! present-day data has no ice. Operates on the aggregated Yelmo-grid smb.
        if (dom%ctl%lim_pd_ice) then
            where(dom%yelmo%dta%pd%H_ice <= 0.0_wp) &
                dom%yelmo%bnd%smb = dom%yelmo%bnd%smb - 4.0_wp
        end if
    end subroutine couple_surface_to_yelmo

    subroutine couple_shelf_to_yelmo(dom)
        ! Basal mass balance + shelf temperature from marine_shelf (grid_shelf ->
        ! Yelmo, conservative).
        type(kryos_domain), intent(inout) :: dom

        real(wp), allocatable :: bmb_y(:,:), Tshlf_y(:,:)
        character(len=256) :: gm, gy

        if (.not. dom%ctl%with_shelf) return

        gm = trim(dom%ctl%grid_shelf)
        gy = trim(dom%ctl%grid_ice)

        call remap(dom, dom%mshlf%now%bmb_shlf, gm, bmb_y,   gy, "con")
        call remap(dom, dom%mshlf%now%T_shlf,   gm, Tshlf_y, gy, "con")
        dom%yelmo%bnd%bmb_shlf = bmb_y
        dom%yelmo%bnd%T_shlf   = Tshlf_y
    end subroutine couple_shelf_to_yelmo

    subroutine step_ice_bnd(dom, ts, coefs)
        ! Prescribed ice ([ice_bnd]): Yelmo's reference ice thickness H_ice_ref
        ! in its fixed and relaxation cells from the ice history at the forcing
        ! time (time_rel), before step_icesheet, which holds the fixed cells at
        ! it and relaxes the others towards it. coefs: the driver's coefficient
        ! registry, advanced to this step (method scaled). No-op unless active.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class),   intent(in)    :: ts
        type(coefs_class),   intent(in)    :: coefs

        call ice_bnd_update(dom%ibnd, dom%yelmo%bnd%H_ice_ref, dom%yelmo%bnd%mask_ice, &
                            dom%yelmo%bnd%tau_relax, ts%time_rel, coefs)
    end subroutine step_ice_bnd

    subroutine step_icesheet(dom, ts)
        ! Advance Yelmo one coupling step on the boundary state assembled by
        ! couple_to_yelmo.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        ! NEGIS: cb_ref from the bed properties + the NEGIS scaling at the
        ! forcing time (Yelmo's till_method = -1, set by domain_negis_init).
        if (dom%ctl%use_negis) &
            call negis_update_cb_ref(dom%yelmo, dom%ngs, ts%time_rel)

        if (.not. dom%ctl%with_ice_sheet) return

        call yelmo_update(dom%yelmo, ts%time)
    end subroutine step_icesheet

    subroutine step_climate(dom, ts, dtt, coefs, init, force)
        ! Run climate on grid_clim, on the dt_clim cadence: the hub geometry
        ! remapped to grid_clim, atmosphere/ocean (and, by backend, surface mass
        ! balance and discharge) produced by the climate backend and its
        ! terms into dom%clim, read by step_surface, step_shelf and
        ! couple_climate_to_yelmo. dtt is the driver's time step. coefs is the
        ! driver-owned coefficient registry (advanced to this step; the climate
        ! sets its derived coefficients in it).
        ! init=.true. marks the cold start: the backend's own initialization
        ! (the REMBO equilibration), with an update regardless of the cadence.
        ! force=.true. updates regardless of the cadence, without that
        ! initialization: the re-derivation on a restart and the refresh after
        ! a cold-start change of the geometry (recon).
        type(kryos_domain),    intent(inout) :: dom
        type(tstep_class),     intent(in)    :: ts
        real(dp),              intent(in)    :: dtt
        type(coefs_class),     intent(inout) :: coefs
        logical,               intent(in), optional :: init
        logical,               intent(in), optional :: force

        real(wp), allocatable :: z_srf_c(:,:), H_ice_c(:,:), z_bed_c(:,:), f_grnd_c(:,:)
        real(wp), allocatable :: z_sl_c(:,:), z_srf_ref_c(:,:)
        character(len=256) :: gc, gh
        character(len=32)  :: mth
        logical :: is_init, is_forced

        if (.not. dom%ctl%with_climate) return

        is_init = .false.
        if (present(init)) is_init = init
        is_forced = .false.
        if (present(force)) is_forced = force

        if (.not. (is_init .or. is_forced .or. cadence_due(ts%time_elapsed, dom%ctl%dt_clim))) return

        gc = trim(dom%ctl%grid_clim)
        gh = trim(dom%ctl%grid_hub)

        ! The hub geometry on grid_clim: averaged over the climate cell when the
        ! hub is finer (conservative), else interpolated.
        mth = remap_method_smooth(dom%ctl%dx_hub, dom%ctl%dx_clim)
        call remap(dom, dom%topo%z_srf,  gh, z_srf_c,  gc, mth)
        call remap(dom, dom%topo%H_ice,  gh, H_ice_c,  gc, mth)
        call remap(dom, dom%topo%z_bed,  gh, z_bed_c,  gc, mth)
        call remap(dom, dom%topo%f_grnd, gh, f_grnd_c, gc, mth)
        call remap(dom, dom%topo%z_sl,   gh, z_sl_c,   gc, mth)
        call remap(dom, dom%topo%ref%z_srf, gh, z_srf_ref_c, gc, mth)

        call climate_update(dom%cl, dom%clim, ts, z_srf_c, H_ice_c, z_bed_c, f_grnd_c, z_sl_c, &
                            z_srf_ref_c, dtt=dtt, mshlf=dom%mshlf, coefs=coefs, init=is_init)
    end subroutine step_climate

    subroutine step_surface(dom, ts, init)
        ! Surface mass balance on grid_surface. Four methods: chion (default; monthly,
        ! needs tas/pr + geometry, with bed and sea level for its land/ice mask),
        ! smbpal (monthly, as chion), smb_simple (needs z_srf + sea-level temperature)
        ! or the climate's own (climate: esm, rembo; nothing to step here). Geometry
        ! comes from the hi-res hub, atmospheric forcing from the climate (grid_clim).
        ! init=.true. runs the snowpack equilibration (chion; smbpal ITM) before the
        ! first update. The result stays on grid_surface in the SMB model's fields
        ! (dom%schn%ann, dom%smb%ann or dom%smbs); couple_surface_to_yelmo lands it
        ! (or the climate's) on the Yelmo grid.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        logical, intent(in), optional    :: init

        real(wp), allocatable :: tas_s(:,:,:), pr_s(:,:,:), z_srf_s(:,:), H_ice_s(:,:)
        real(wp), allocatable :: tsl_s(:,:), Href_s(:,:), z_bed_s(:,:), z_sl_s(:,:)
        character(len=256) :: gc, gs, gh, gy
        character(len=32)  :: mth
        logical :: is_init

        if (.not. dom%ctl%with_surface) return

        is_init = .false.
        if (present(init)) is_init = init

        gc = trim(dom%ctl%grid_clim)
        gs = trim(dom%ctl%grid_surface)
        gh = trim(dom%ctl%grid_hub)
        gy = trim(dom%ctl%grid_ice)

        ! The hub geometry on grid_surface: averaged over the surface cell when the
        ! hub is finer (conservative), else interpolated.
        mth = remap_method_smooth(dom%ctl%dx_hub, dom%ctl%dx_surface)

        ! surface_method = "climate": nothing to do, the climate has the smb
        ! (couple_surface_to_yelmo lands it).
        if (trim(dom%ctl%surface_method) == "smb_simple") then
            ! smb_simple: surface elevation + sea-level temperature, masked to the
            ! reference ice extent (refreshed each call in case H_ice_ref changed).
            call remap(dom, dom%topo%z_srf,          gh, z_srf_s, gs, mth)
            call remap(dom, dom%clim%now%tsl_ann,     gc, tsl_s,   gs, "bilin")
            call remap(dom, dom%yelmo%bnd%H_ice_ref,  gy, Href_s,  gs, "bilin")
            call smb_simple_set_mask(dom%smbs, Href_s)
            call smb_simple_update(dom%smbs, z_srf_s, tsl_s)
        else if (trim(dom%ctl%surface_method) == "chion") then
            ! chion (monthly -> daily annual cycle)
            call remap(dom, dom%topo%z_srf,   gh, z_srf_s, gs, mth)
            call remap(dom, dom%topo%H_ice,   gh, H_ice_s, gs, mth)
            call remap(dom, dom%topo%z_bed,   gh, z_bed_s, gs, mth)
            call remap(dom, dom%topo%z_sl,    gh, z_sl_s,  gs, mth)
            call remap_atm_to_surface(dom, z_srf_s, tas_s, pr_s)
            if (is_init) call surface_chion_spinup(dom%schn, tas_s, pr_s, z_srf_s, z_bed_s, H_ice_s, &
                                                   z_sl_s, ts%time_rel)
            call surface_chion_update(dom%schn, tas_s, pr_s, z_srf_s, z_bed_s, H_ice_s, z_sl_s, &
                                      ts%time_rel)
        else if (trim(dom%ctl%surface_method) == "smbpal") then
            ! smbpal (monthly)
            call remap(dom, dom%topo%z_srf,  gh, z_srf_s, gs, mth)
            call remap(dom, dom%topo%H_ice,  gh, H_ice_s, gs, mth)
            call remap_atm_to_surface(dom, z_srf_s, tas_s, pr_s)
            if (is_init .and. trim(dom%smb%par%abl_method) == "itm") then
                call smbpal_update_monthly_equil(dom%smb, tas_s, pr_s, z_srf_s, H_ice_s, ts%time_rel)
            end if
            call smbpal_update_monthly(dom%smb, tas_s, pr_s, z_srf_s, H_ice_s, ts%time_rel)
        end if
    end subroutine step_surface

    subroutine remap_atm_to_surface(dom, z_srf_s, tas_s, pr_s)
        ! Monthly tas and pr from grid_clim onto grid_surface at its surface
        ! elevation z_srf_s: interpolated (bilinear), then moved from the surface
        ! they were computed for (clim z_srf, interpolated) to z_srf_s with the
        ! backend's monthly lapse rate and precipitation sensitivity:
        !   tas_s = tas - lapse*dz,  pr_s = pr*exp(-beta_p*lapse*dz),  dz = z_srf_s - z_srf.
        ! For a bilinear remap this equals interpolating the sea-level fields and
        ! lifting them to z_srf_s. On the same grid dz = 0: a plain copy.
        type(kryos_domain),    intent(inout) :: dom
        real(wp),              intent(in)    :: z_srf_s(:,:)
        real(wp), allocatable, intent(inout) :: tas_s(:,:,:), pr_s(:,:,:)

        real(wp), allocatable :: zc_s(:,:), bp_s(:,:), dz(:,:)
        character(len=256) :: gc, gs
        integer :: m

        gc = trim(dom%ctl%grid_clim)
        gs = trim(dom%ctl%grid_surface)

        call remap(dom, dom%clim%now%tas,    gc, tas_s, gs, "bilin")
        call remap(dom, dom%clim%now%pr,     gc, pr_s,  gs, "bilin")
        if (trim(gc) == trim(gs)) return

        call remap(dom, dom%clim%now%z_srf,  gc, zc_s,  gs, "bilin")
        call remap(dom, dom%clim%now%beta_p, gc, bp_s,  gs, "bilin")
        dz = z_srf_s - zc_s
        do m = 1, size(tas_s,3)
            tas_s(:,:,m) = tas_s(:,:,m) - dom%clim%lapse(m)*dz
            pr_s(:,:,m)  = pr_s(:,:,m) * exp(-bp_s*dom%clim%lapse(m)*dz)
        end do

    end subroutine remap_atm_to_surface

    subroutine couple_yelmo_to_htopo(dom)
        ! Land Yelmo's state on the hub's current geometry. On the Yelmo grid
        ! there is no finer information, so the hub mirrors Yelmo (including its
        ! fractional grounding). On a finer hub, the hub keeps its hi-res reference
        ! and adds Yelmo's anomalies, refined bilinearly: the bed displacement
        ! (z_bed - z_bed_ref) and the change in ice thickness from the hub
        ! reference as Yelmo received it (conservative); htopo_update recomputes
        ! the grounding and the surface on the hub. Static masks are not refreshed.
        type(kryos_domain), intent(inout) :: dom

        real(wp), allocatable :: H_ice_ref_y(:,:), dz_bed_h(:,:), dH_ice_h(:,:), z_sl_h(:,:)
        character(len=256) :: gh, gy

        gh = trim(dom%ctl%grid_hub)
        gy = trim(dom%ctl%grid_ice)

        if (trim(gh) == trim(gy)) then
            dom%topo%H_ice  = dom%yelmo%tpo%now%H_ice
            dom%topo%z_bed  = dom%yelmo%bnd%z_bed
            dom%topo%f_grnd = dom%yelmo%tpo%now%f_grnd
            dom%topo%z_sl   = dom%yelmo%bnd%z_sl
            dom%topo%z_srf  = dom%yelmo%tpo%now%z_srf
        else
            call remap(dom, dom%topo%ref%H_ice, gh, H_ice_ref_y, gy, "con")
            call remap(dom, dom%yelmo%bnd%z_bed - dom%yelmo%bnd%z_bed_ref, gy, dz_bed_h, gh, "bilin")
            call remap(dom, dom%yelmo%tpo%now%H_ice - H_ice_ref_y,         gy, dH_ice_h, gh, "bilin")
            call remap(dom, dom%yelmo%bnd%z_sl,                            gy, z_sl_h,   gh, "bilin")
            call htopo_update(dom%topo, dz_bed_h, dH_ice_h, z_sl_h)
        end if
    end subroutine couple_yelmo_to_htopo

    subroutine step_shelf(dom, ts)
        ! Run marine_shelf on its own grid: geometry/masks from the hub, ocean
        ! forcing from the climate, as depth profiles (interpolated to the shelf
        ! base here) or already at the shelf base (has_ocn_shelf). The outputs stay on grid_shelf (in dom%mshlf%now);
        ! couple_shelf_to_yelmo lands bmb_shlf / T_shlf on the Yelmo grid.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        real(wp), allocatable :: H_ice_m(:,:), z_bed_m(:,:), f_grnd_m(:,:), z_sl_m(:,:)
        real(wp), allocatable :: z_srf_m(:,:)
        real(wp), allocatable :: to_m(:,:,:), so_m(:,:,:), dto_m(:,:,:), dto_y(:,:,:)
        character(len=256) :: gm, gh, gc

        if (.not. dom%ctl%with_shelf) return

        gm = trim(dom%ctl%grid_shelf)
        gh = trim(dom%ctl%grid_hub)
        gc = trim(dom%ctl%grid_clim)

        ! geometry: hub -> mshlf grid
        call remap(dom, dom%topo%H_ice,   gh, H_ice_m,   gm, "bilin")
        call remap(dom, dom%topo%z_bed,   gh, z_bed_m,   gm, "bilin")
        call remap(dom, dom%topo%f_grnd,  gh, f_grnd_m,  gm, "bilin")
        call remap(dom, dom%topo%z_sl,    gh, z_sl_m,    gm, "bilin")
        call remap(dom, dom%topo%z_srf,   gh, z_srf_m,   gm, "bilin")

        if (dom%clim%has_ocn_shelf) then
            ! Ocean already at the shelf base (esm): grid_clim -> mshlf grid.
            call remap(dom, dom%clim%now%T_shlf,  gc, dom%mshlf%now%T_shlf,  gm, "bilin")
            call remap(dom, dom%clim%now%S_shlf,  gc, dom%mshlf%now%S_shlf,  gm, "bilin")
            call remap(dom, dom%clim%now%dT_shlf, gc, dom%mshlf%now%dT_shlf, gm, "bilin")
            call remap(dom, dom%clim%now%dS_shlf, gc, dom%mshlf%now%dS_shlf, gm, "bilin")
            call marshelf_update(dom%mshlf, H_ice_m, z_bed_m, f_grnd_m, &
                    z_sl_m, dx=dom%ctl%dx_shelf, z_srf=z_srf_m)
            return
        end if

        ! ocean forcing (3D): depth profiles (grid_clim) -> mshlf grid
        call remap(dom, dom%clim%now%to_ann, gc, to_m, gm, "bilin")
        call remap(dom, dom%clim%now%so_ann, gc, so_m, gm, "bilin")
        dto_y = dom%clim%now%to_ann - dom%clim%ref%to_ann
        call remap(dom, dto_y, gc, dto_m, gm, "bilin")

        ! run marine_shelf on grid_shelf
        call marshelf_update_shelf(dom%mshlf, H_ice_m, z_bed_m, f_grnd_m, z_sl_m, &
                dom%ctl%dx_shelf, dom%clim%now%depth, to_m, so_m, dto_ann=dto_m)
        call marshelf_update(dom%mshlf, H_ice_m, z_bed_m, f_grnd_m, &
                z_sl_m, dx=dom%ctl%dx_shelf, z_srf=z_srf_m)
    end subroutine step_shelf

end module kryos_coupling
