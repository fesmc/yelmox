module kryos_coupling
    ! The per-step coupling primitives of a kryos_domain. Each step_* advances
    ! one component on its own grid, reading geometry from the hi-res hub; each
    ! couple_*_to_yelmo lands one component's output on the Yelmo grid, before
    ! Yelmo runs. Every primitive is a no-op when its component is inactive.

    use timestepping, only : tstep_class
    use yelmo,        only : wp, yelmo_update
    use marine_shelf, only : marshelf_update, marshelf_update_shelf
    use fastisostasy, only : isos_update, bsl_class
    use yelmox_climate, only : climate_update
    use smbpal,       only : smbpal_update_monthly, smbpal_update_monthly_equil
    use smb_simple_m, only : smb_simple_set_mask, smb_simple_update
    use htopo,        only : htopo_update
    use ice_optimization, only : optimize_set_transient_param, optimize_cb_ref, optimize_tf_corr
    use kryos,        only : kryos_domain, remap, remap_method_smooth, cadence_due
    use kryos_regions, only : negis_update_cb_ref, calc_glacial_smb
    use kryos_forcing, only : tsforcing_class

    implicit none
    private

    public :: step_spinup_tuning, step_isostasy, step_icesheet, step_climate, step_marine_shelf
    public :: step_smb, refresh_hub, update_climate
    public :: couple_to_yelmo
    public :: couple_isostasy_to_yelmo, couple_smb_to_yelmo, couple_marine_to_yelmo
    public :: couple_climate_to_yelmo
    public :: check_isostasy_reference

contains

    subroutine step_spinup_tuning(dom, ts)
        ! Spin-up tuning (equil_method == "opt"): ramp the topography relaxation
        ! timescale, then nudge the basal-friction field cb_ref and the marine
        ! thermal-forcing correction tf_corr toward present-day observations.
        !
        ! cb_ref is a Yelmo-grid control, optimized in place. tf_corr lives on the
        ! marine_shelf grid; the observational targets (H_ice/H_grnd) live on the
        ! Yelmo grid, so the correction is lifted to the Yelmo grid (tf_corr_y),
        ! optimized there, and remapped back to the shelf grid. At identity grids
        ! both remaps are copies.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        real(wp), allocatable :: tf_corr_y(:,:), tf_corr_m(:,:)
        character(len=256) :: gm, gy

        if (trim(dom%ctl%equil_method) /= "opt") return

        gm = trim(dom%ctl%grid_mshlf)
        gy = trim(dom%ctl%grid_ice)

        ! Topography relaxation ramp (gl + grounding-zone relaxing while active).
        if (ts%time_elapsed <= dom%opt%rel_time2) then
            call optimize_set_transient_param(dom%opt%rel_tau, ts%time_elapsed, &
                    time1=dom%opt%rel_time1, time2=dom%opt%rel_time2, &
                    p1=dom%opt%rel_tau1, p2=dom%opt%rel_tau2, m=dom%opt%rel_m)
            dom%yelmo%tpo%par%topo_rel_tau = dom%opt%rel_tau
            dom%yelmo%tpo%par%topo_rel     = 4
        else
            dom%yelmo%tpo%par%topo_rel = 0
        end if

        ! Basal friction (cb_ref) optimization -- Yelmo grid, in place.
        if (dom%opt%opt_cf .and. ts%time_elapsed >= dom%opt%cf_time_init &
                            .and. ts%time_elapsed <= dom%opt%cf_time_end) then
            call optimize_cb_ref(dom%yelmo%dyn%now%cb_ref, dom%yelmo%tpo%now%H_ice, &
                    dom%yelmo%tpo%now%dHidt, dom%yelmo%bnd%z_bed, dom%yelmo%bnd%z_sl, &
                    dom%yelmo%dyn%now%ux_s, dom%yelmo%dyn%now%uy_s, &
                    dom%yelmo%dta%pd%H_ice, dom%yelmo%dta%pd%uxy_s, dom%yelmo%dta%pd%H_grnd, &
                    dom%opt%cf_min, dom%opt%cf_max, dom%yelmo%tpo%par%dx, &
                    dom%opt%sigma_err, dom%opt%sigma_vel, dom%opt%tau_c, dom%opt%H0, &
                    dt=dom%ctl%dtt, fill_method=dom%opt%fill_method, fill_dist=dom%opt%sigma_err, &
                    cb_tgt=dom%yelmo%dyn%now%cb_tgt)
        end if

        ! Thermal-forcing correction (tf_corr) optimization -- lift shelf-grid
        ! correction to the Yelmo grid, optimize against Yelmo-grid targets, remap
        ! back. tf_corr persists on the shelf grid (in mshlf, incl. its restart).
        if (dom%opt%opt_tf .and. ts%time_elapsed >= dom%opt%tf_time_init &
                            .and. ts%time_elapsed <= dom%opt%tf_time_end) then
            call remap(dom, dom%mshlf%now%tf_corr, gm, tf_corr_y, gy, "con")
            call optimize_tf_corr(tf_corr_y, dom%yelmo%tpo%now%H_ice, dom%yelmo%tpo%now%H_grnd, &
                    dom%yelmo%tpo%now%dHidt, dom%yelmo%dta%pd%H_ice, dom%yelmo%dta%pd%H_grnd, &
                    dom%opt%H_grnd_lim, dom%yelmo%bnd%basins, dom%opt%basin_fill, &
                    dom%opt%tau_m, dom%opt%m_temp, dom%opt%tf_min, dom%opt%tf_max, &
                    dom%yelmo%tpo%par%dx, sigma=dom%opt%tf_sigma, dt=dom%ctl%dtt)
            call remap(dom, tf_corr_y, gy, tf_corr_m, gm, "bilin")
            dom%mshlf%now%tf_corr = tf_corr_m
        end if
    end subroutine step_spinup_tuning

    subroutine step_isostasy(dom, ts, bsl)
        ! Run isostasy on its own grid: ice load from Yelmo (bilin). The bedrock /
        ! sea-surface outputs stay on grid_isos (in dom%isos%out); they are landed
        ! on the Yelmo grid by couple_isostasy_to_yelmo (in step_icesheet, before
        ! yelmo_update). Assumes grid_isos is at least as fine as grid_ice
        ! (identity when equal). bsl is the shared, driver-owned sea level (already
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
        call couple_smb_to_yelmo(dom)
        call couple_marine_to_yelmo(dom)
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
            stop
        end if
    end subroutine check_isostasy_reference

    subroutine couple_smb_to_yelmo(dom)
        ! Surface mass balance + surface temperature from the active SMB model
        ! (grid_smb -> Yelmo, conservative), with the we->ie unit scaling and the
        ! optional Greenland modifications. The producing step (step_smb,
        ! or a flavor climate step) leaves smb/tsrf on grid_smb in the SMB model's
        ! own fields; this coupler is the single place that lands them on Yelmo.
        type(kryos_domain), intent(inout) :: dom

        real(wp), allocatable :: smb_y(:,:), tsrf_y(:,:), ta_y(:,:), ta_pd_y(:,:)
        character(len=256) :: gs, gc, gy

        if (.not. dom%ctl%with_climate) return

        gs = trim(dom%ctl%grid_smb)
        gc = trim(dom%ctl%grid_clim)
        gy = trim(dom%ctl%grid_ice)

        if (trim(dom%ctl%smb_method) == "smb_simple") then
            call remap(dom, dom%smbs%smb,   gs, smb_y,  gy, "con")
            call remap(dom, dom%smbs%t_srf, gs, tsrf_y, gy, "con")
        else
            call remap(dom, dom%smb%ann%smb,  gs, smb_y,  gy, "con")
            call remap(dom, dom%smb%ann%tsrf, gs, tsrf_y, gy, "con")
        end if

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
    end subroutine couple_smb_to_yelmo

    subroutine couple_marine_to_yelmo(dom)
        ! Basal mass balance + shelf temperature from marine_shelf (grid_mshlf ->
        ! Yelmo, conservative).
        type(kryos_domain), intent(inout) :: dom

        real(wp), allocatable :: bmb_y(:,:), Tshlf_y(:,:)
        character(len=256) :: gm, gy

        if (.not. dom%ctl%with_marine_shelf) return

        gm = trim(dom%ctl%grid_mshlf)
        gy = trim(dom%ctl%grid_ice)

        call remap(dom, dom%mshlf%now%bmb_shlf, gm, bmb_y,   gy, "con")
        call remap(dom, dom%mshlf%now%T_shlf,   gm, Tshlf_y, gy, "con")
        dom%yelmo%bnd%bmb_shlf = bmb_y
        dom%yelmo%bnd%T_shlf   = Tshlf_y
    end subroutine couple_marine_to_yelmo

    subroutine step_icesheet(dom, ts)
        ! Advance Yelmo one coupling step on the boundary state assembled by
        ! couple_to_yelmo.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        ! NEGIS: update cb_ref from bed properties + NEGIS scaling.
        if (dom%ngs%use_negis_par) &
            call negis_update_cb_ref(dom%yelmo, dom%ngs, ts%time)

        if (.not. dom%ctl%with_ice_sheet) return
        if (ts%n == 0 .and. dom%yelmo%par%use_restart) return

        call yelmo_update(dom%yelmo, ts%time)
    end subroutine step_icesheet

    subroutine step_climate(dom, ts, tsf)
        ! Run climate on grid_clim, on the dt_clim cadence: geometry from the hub,
        ! atmosphere/ocean (and, by backend, surface mass balance and discharge)
        ! produced by the climate backend into dom%clim, read by step_smb,
        ! step_marine_shelf and couple_climate_to_yelmo. tsf (optional) is the
        ! driver-owned transient forcing; see update_climate.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        type(tsforcing_class), intent(in), optional :: tsf

        if (.not. dom%ctl%with_climate) return

        if (cadence_due(ts%time_elapsed, dom%ctl%dt_clim)) call update_climate(dom, ts, tsf=tsf)
    end subroutine step_climate

    subroutine update_climate(dom, ts, tsf, init)
        ! One climate-backend update on grid_clim: the hub geometry remapped to
        ! grid_clim, with the transient forcing (tsf) when given; the backend
        ! applies it in its own way. init marks the cold start.
        type(kryos_domain),    intent(inout) :: dom
        type(tstep_class),     intent(in)    :: ts
        type(tsforcing_class), intent(in), optional :: tsf
        logical,               intent(in), optional :: init

        real(wp), allocatable :: z_srf_c(:,:), H_ice_c(:,:), z_bed_c(:,:), f_grnd_c(:,:)
        real(wp), allocatable :: z_sl_c(:,:), z_srf_ref_c(:,:), basins_c(:,:)
        character(len=256) :: gc, gh

        gc = trim(dom%ctl%grid_clim)
        gh = trim(dom%ctl%grid_hub)

        call remap(dom, dom%topo%z_srf,  gh, z_srf_c,  gc, "bilin")
        call remap(dom, dom%topo%H_ice,  gh, H_ice_c,  gc, "bilin")
        call remap(dom, dom%topo%z_bed,  gh, z_bed_c,  gc, "bilin")
        call remap(dom, dom%topo%f_grnd, gh, f_grnd_c, gc, "bilin")
        call remap(dom, dom%topo%z_sl,   gh, z_sl_c,   gc, "bilin")
        call remap(dom, dom%topo%z_srf_ref, gh, z_srf_ref_c, gc, "bilin")
        call remap(dom, dom%topo%basins, gh, basins_c, gc, "nn")

        call climate_update(dom%cl, dom%clim, ts, z_srf_c, H_ice_c, z_bed_c, f_grnd_c, z_sl_c, &
                            z_srf_ref_c, basins_c, domain=dom%ctl%domain, dx=dom%ctl%dx_clim, &
                            dtt=dom%ctl%dtt, mshlf=dom%mshlf, tsf=tsf, init=init)
    end subroutine update_climate

    subroutine step_smb(dom, ts, init)
        ! Surface mass balance on grid_smb. Two methods: smbpal (default; monthly,
        ! needs tas/pr + geometry) or smb_simple (needs z_srf + sea-level
        ! temperature). Geometry comes from the hi-res hub, atmospheric forcing from
        ! snapclim (grid_clim). init=.true. runs the smbpal ITM equilibration before
        ! the first update. The result stays on grid_smb in the SMB model's fields
        ! (dom%smb%ann or dom%smbs); couple_smb_to_yelmo lands it on the Yelmo grid.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        logical, intent(in), optional    :: init

        real(wp), allocatable :: tas_s(:,:,:), pr_s(:,:,:), z_srf_s(:,:), H_ice_s(:,:)
        real(wp), allocatable :: tsl_s(:,:), Href_s(:,:)
        real(wp), allocatable :: smb_s(:,:), tsrf_s(:,:)
        character(len=256) :: gc, gs, gh, gy
        logical :: is_init

        if (.not. dom%ctl%with_climate) return

        is_init = .false.
        if (present(init)) is_init = init

        gc = trim(dom%ctl%grid_clim)
        gs = trim(dom%ctl%grid_smb)
        gh = trim(dom%ctl%grid_hub)
        gy = trim(dom%ctl%grid_ice)

        if (trim(dom%ctl%smb_method) == "climate") then
            ! The climate's own surface mass balance and surface temperature, at
            ! the current surface.
            call remap(dom, dom%clim%now%smb,  gc, smb_s,  gs, "bilin")
            call remap(dom, dom%clim%now%tsrf, gc, tsrf_s, gs, "bilin")
            dom%smb%ann%smb  = smb_s
            dom%smb%ann%tsrf = tsrf_s
        else if (trim(dom%ctl%smb_method) == "smb_simple") then
            ! smb_simple: surface elevation + sea-level temperature, masked to the
            ! reference ice extent (refreshed each call in case H_ice_ref changed).
            call remap(dom, dom%topo%z_srf,          gh, z_srf_s, gs, "bilin")
            call remap(dom, dom%clim%now%tsl_ann,     gc, tsl_s,   gs, "bilin")
            call remap(dom, dom%yelmo%bnd%H_ice_ref,  gy, Href_s,  gs, "bilin")
            call smb_simple_set_mask(dom%smbs, Href_s)
            call smb_simple_update(dom%smbs, z_srf_s, tsl_s)
        else
            ! smbpal (monthly)
            call remap(dom, dom%clim%now%tas, gc, tas_s, gs, "bilin")
            call remap(dom, dom%clim%now%pr,  gc, pr_s,  gs, "bilin")
            call remap(dom, dom%topo%z_srf,  gh, z_srf_s, gs, "bilin")
            call remap(dom, dom%topo%H_ice,  gh, H_ice_s, gs, "bilin")
            if (is_init .and. trim(dom%smb%par%abl_method) == "itm") then
                call smbpal_update_monthly_equil(dom%smb, tas_s, pr_s, z_srf_s, H_ice_s, &
                        ts%time_rel, time_equil=100.0_wp)
            end if
            call smbpal_update_monthly(dom%smb, tas_s, pr_s, z_srf_s, H_ice_s, ts%time_rel)
        end if
    end subroutine step_smb

    subroutine refresh_hub(dom)
        ! Refresh the hub's current geometry from the models. On the Yelmo grid
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
            call remap(dom, dom%topo%H_ice_ref, gh, H_ice_ref_y, gy, "con")
            call remap(dom, dom%yelmo%bnd%z_bed - dom%yelmo%bnd%z_bed_ref, gy, dz_bed_h, gh, "bilin")
            call remap(dom, dom%yelmo%tpo%now%H_ice - H_ice_ref_y,         gy, dH_ice_h, gh, "bilin")
            call remap(dom, dom%yelmo%bnd%z_sl,                            gy, z_sl_h,   gh, "bilin")
            call htopo_update(dom%topo, dz_bed_h, dH_ice_h, z_sl_h)
        end if
    end subroutine refresh_hub

    subroutine step_marine_shelf(dom, ts)
        ! Run marine_shelf on its own grid: geometry/masks from the hub, ocean
        ! forcing from the climate, as depth profiles (interpolated to the shelf
        ! base here) or already at the shelf base (has_ocn_shelf). The outputs stay on grid_mshlf (in dom%mshlf%now);
        ! couple_marine_to_yelmo lands bmb_shlf / T_shlf on the Yelmo grid.
        type(kryos_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        real(wp), allocatable :: H_ice_m(:,:), z_bed_m(:,:), f_grnd_m(:,:), z_sl_m(:,:)
        real(wp), allocatable :: z_srf_m(:,:)
        real(wp), allocatable :: regions_m(:,:), basins_m(:,:)
        real(wp), allocatable :: to_m(:,:,:), so_m(:,:,:), dto_m(:,:,:), dto_y(:,:,:)
        character(len=256) :: gm, gh, gc

        if (.not. dom%ctl%with_marine_shelf) return

        gm = trim(dom%ctl%grid_mshlf)
        gh = trim(dom%ctl%grid_hub)
        gc = trim(dom%ctl%grid_clim)

        ! geometry + masks: hub -> mshlf grid
        call remap(dom, dom%topo%H_ice,   gh, H_ice_m,   gm, "bilin")
        call remap(dom, dom%topo%z_bed,   gh, z_bed_m,   gm, "bilin")
        call remap(dom, dom%topo%f_grnd,  gh, f_grnd_m,  gm, "bilin")
        call remap(dom, dom%topo%z_sl,    gh, z_sl_m,    gm, "bilin")
        call remap(dom, dom%topo%z_srf,   gh, z_srf_m,   gm, "bilin")
        call remap(dom, dom%topo%regions, gh, regions_m, gm, "nn")
        call remap(dom, dom%topo%basins,  gh, basins_m,  gm, "nn")

        if (dom%clim%has_ocn_shelf) then
            ! Ocean already at the shelf base (esm): grid_clim -> mshlf grid.
            call remap(dom, dom%clim%now%T_shlf,  gc, dom%mshlf%now%T_shlf,  gm, "bilin")
            call remap(dom, dom%clim%now%S_shlf,  gc, dom%mshlf%now%S_shlf,  gm, "bilin")
            call remap(dom, dom%clim%now%dT_shlf, gc, dom%mshlf%now%dT_shlf, gm, "bilin")
            call remap(dom, dom%clim%now%dS_shlf, gc, dom%mshlf%now%dS_shlf, gm, "bilin")
            call marshelf_update(dom%mshlf, H_ice_m, z_bed_m, f_grnd_m, regions_m, basins_m, &
                    z_sl_m, dx=dom%ctl%dx_mshlf, z_srf=z_srf_m)
            return
        end if

        ! ocean forcing (3D): depth profiles (grid_clim) -> mshlf grid
        call remap(dom, dom%clim%now%to_ann, gc, to_m, gm, "bilin")
        call remap(dom, dom%clim%now%so_ann, gc, so_m, gm, "bilin")
        dto_y = dom%clim%now%to_ann - dom%clim%ref%to_ann
        call remap(dom, dto_y, gc, dto_m, gm, "bilin")

        ! run marine_shelf on grid_mshlf
        call marshelf_update_shelf(dom%mshlf, H_ice_m, z_bed_m, f_grnd_m, basins_m, z_sl_m, &
                dom%ctl%dx_mshlf, dom%clim%now%depth, to_m, so_m, dto_ann=dto_m)
        call marshelf_update(dom%mshlf, H_ice_m, z_bed_m, f_grnd_m, regions_m, basins_m, &
                z_sl_m, dx=dom%ctl%dx_mshlf, z_srf=z_srf_m)
    end subroutine step_marine_shelf

end module kryos_coupling
