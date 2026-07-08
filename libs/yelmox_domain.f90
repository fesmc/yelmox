module yelmox_domain
    ! Multigrid yelmox: one region's full model state bundled as an ice_domain,
    ! advanced by composable step_* primitives.
    !
    ! The hi-res topography hub (dom%topo) is the geometry source of truth: it is
    ! refreshed from the prognostic models each step (refresh_htopo), and the
    ! coupling steps remap fields to/from it through the domain's coupler
    ! (dom%cpl). Each helper model runs on its own grid, named in domain_ctl; a
    ! component's grid is just a string, so e.g. marine_shelf can run on the
    ! hi-res hub grid or on the Yelmo grid simply by changing grid_mshlf.
    !
    ! Status: isostasy, climate, smb and marine_shelf are each coupled through the
    ! hub on their own configurable grids (grid_isos/grid_clim/grid_smb/grid_mshlf).

    use nml,          only : nml_read
    use ncio
    use timestepping, only : tstep_class, tstep_init
    use coords,       only : grid_class, grid_cdo_read_desc
    use yelmo,        only : yelmo_class, wp, yelmo_init, yelmo_update, yelmo_update_equil, &
                             yelmo_init_state, yelmo_init_topo, yelmo_print_bound, &
                             yelmo_restart_write, yelmo_restart_read, &
                             yelmo_regions_init, yelmo_region_init, yelmo_regions_update, &
                             yelmo_write_init, yelmo_write_step, yelmo_regions_write
    use yelmo_defs,   only : MASK_ICE_NONE, MASK_ICE_FIXED, MASK_ICE_DYNAMIC
    use yelmo_tools,  only : smooth_gauss_2D
    use basal_dragging, only : calc_cb_ref
    use ice_sub_regions, only : get_ice_sub_region
    use yelmo_topography, only : calc_ytopo_diagnostic
    use yelmo_io,         only : yelmo_restart_read_topo_bnd
    use marine_shelf, only : marshelf_class, marshelf_init, marshelf_update, &
                             marshelf_update_shelf, marshelf_restart_write, &
                             marshelf_restart_read
    use fastisostasy, only : isos_class, isos_init, isos_update, isos_init_ref, &
                             isos_init_state, isos_restart_write, isos_restart_read, &
                             bsl_class, bsl_update, bsl_restart_read, bsl_restart_write
    use snapclim,     only : snapclim_class, snapclim_init, snapclim_update
    use smbpal,       only : smbpal_class, smbpal_init, smbpal_update_monthly, &
                             smbpal_update_monthly_equil, smbpal_restart_write, smbpal_restart_read
    use smb_simple_m, only : smb_simple_class, smb_simple_init, smb_simple_set_mask, &
                             smb_simple_update
    use ice_optimization, only : ice_opt_params, optimize_par_load, &
                             optimize_set_transient_param, optimize_cb_ref, optimize_tf_corr
    use sediments,    only : sediments_class, sediments_init
    use geothermal,   only : geothermal_class, geothermal_init
    use htopo,        only : htopo_class, htopo_init, htopo_write_init, htopo_write_step
    use coupler,      only : coupler_class, coupler_init, coupler_prime, cpl_remap => remap

    implicit none
    private

    character(len=*), parameter :: MAP_FLDR = "maps"

    type domain_ctl
        ! Parameter file (kept for sub-steps that reload from it, e.g. LGM startup).
        character(len=512) :: path_par = ""
        ! Shared timeline (read from the driver's timeline group, e.g. [ctrl]).
        character(len=56) :: tstep_method = "const"
        real(wp) :: dtt         = 10.0_wp
        ! Cadences + methods ([coupling]).
        real(wp) :: dt_restart  = 0.0_wp
        real(wp) :: dt_clim     = 10.0_wp   ! [yr] snapclim snapshot update frequency
        character(len=56) :: equil_method = "none"
        character(len=56) :: smb_method   = "smbpal"

        ! Domain-specific startup / physics switches ([coupling]; Greenland only).
        logical :: greenland_init_marine_H = .false.   ! impose LGM-like marine ice at start
        logical :: scale_glacial_smb       = .false.   ! reduce negative glacial smb (Greenland)
        logical :: lim_pd_ice              = .false.   ! extra melt outside PD ice extent (Greenland/rembo)
        logical :: use_negis               = .false.   ! NEGIS cb_ref modification (Greenland)

        ! Which components are active in this domain's coupling sequence.
        logical :: with_ice_sheet    = .true.
        logical :: with_isostasy     = .true.
        logical :: with_marine_shelf = .true.
        logical :: with_climate      = .true.

        ! Domain + grid names (the source of truth for remap keys).
        character(len=256) :: domain     = ""   ! e.g. "Antarctica"
        character(len=256) :: grid_name  = ""   ! hi-res reference (htopo), highest res
        character(len=256) :: grid_yelmo = ""   ! Yelmo grid
        character(len=256) :: grid_mshlf = ""   ! marine-shelf grid ([coupling]; default = grid_name)
        real(wp) :: dx_mshlf = 0.0_wp           ! marine-shelf grid spacing (Yelmo dx units)
        character(len=256) :: grid_isos = ""    ! isostasy grid ([coupling]; default = grid_yelmo)
        real(wp) :: dx_isos = 0.0_wp            ! isostasy grid spacing in x (Yelmo dx units)
        real(wp) :: dy_isos = 0.0_wp            ! isostasy grid spacing in y (Yelmo dy units)
        character(len=256) :: grid_clim = ""    ! climate grid ([coupling]; default = grid_yelmo)
        real(wp) :: dx_clim = 0.0_wp            ! climate grid spacing (Yelmo dx units)
        character(len=256) :: grid_smb = ""     ! smb grid ([coupling]; default = grid_clim)

        ! Restart bundle folder ([coupling]); "None" = cold start.
        character(len=512) :: restart = "None"

        ! Per-module output switches ([output]); each module -> its own file.
        logical :: write_yelmo = .true.
        logical :: write_isos  = .true.
        logical :: write_mshlf = .true.
        logical :: write_smb   = .true.
        logical :: write_snap  = .true.
        logical :: write_htopo = .true.
    end type domain_ctl

    ! NEGIS (Northeast Greenland Ice Stream) cb_ref modification parameters.
    type negis_params
        logical  :: use_negis_par = .false.
        real(wp) :: cf_0    = 1.0_wp
        real(wp) :: cf_1    = 1.0_wp
        real(wp) :: cf_centre = 1.0_wp
        real(wp) :: cf_north  = 1.0_wp
        real(wp) :: cf_south  = 1.0_wp
        real(wp) :: cf_x    = 1.0_wp
    end type negis_params

    type ice_domain
        type(yelmo_class)      :: yelmo
        type(marshelf_class)   :: mshlf
        type(isos_class)       :: isos
        type(snapclim_class)   :: snp
        type(smbpal_class)     :: smb
        type(smb_simple_class) :: smbs    ! alternative SMB (smb_method="smb_simple")
        type(sediments_class)  :: sed
        type(geothermal_class) :: gthrm
        type(htopo_class)      :: topo    ! hi-res geometry reference hub
        type(coupler_class)    :: cpl     ! this region's grid resolution + map cache
        type(ice_opt_params)   :: opt     ! basal-friction / thermal-forcing optimization
        type(negis_params)     :: ngs     ! Greenland NEGIS cb_ref modification
        type(domain_ctl)       :: ctl
    end type ice_domain

    public :: domain_ctl, ice_domain
    public :: timeline_init, tstep_due
    public :: domain_init, domain_regions_init, domain_init_state, yelmox_step
    public :: domain_startup, bsl_startup, run_restart_write
    public :: domain_restart_write, domain_restart_read, restart_bundle_dir, restart_bundle_mkdir
    public :: domain_write_init, domain_write_step, domain_write_1D
    public :: step_isostasy, step_icesheet, step_climate, refresh_htopo, step_marine_shelf
    public :: step_optimize, domain_update_smb
    public :: couple_isostasy_to_yelmo, couple_smb_to_yelmo, couple_marine_to_yelmo
    ! Exposed so flavor-specific drivers (e.g. the ESM driver) can remap between
    ! the hub/Yelmo grids and their own forcing grid via the domain's coupler.
    public :: remap

    ! Domain-level remap: identity-copy when src == dst, else remap via the coupler.
    interface remap
        module procedure remap_2D, remap_3D
    end interface remap

contains

    subroutine timeline_init(ts, dtt, path_par, group, time_ref, cal)
        ! Read the run's shared timeline ([<group>]: tstep_method, tstep_const,
        ! time_init, time_end, dtt) and initialize the driver-owned timestepper.
        ! The same group name is passed to domain_init (timeline_group), so the
        ! domain reads the timeline values it needs itself -- drivers never
        ! inject them. time_ref sets the calendar reference (default 1950.0);
        ! cal=.true. applies tstep_const as a calendar constant (const_cal, the
        ! ESM convention) instead of a relative one (const_rel).
        type(tstep_class), intent(out) :: ts
        real(wp),          intent(out) :: dtt
        character(len=*),  intent(in)  :: path_par
        character(len=*),  intent(in)  :: group
        real(wp), optional, intent(in) :: time_ref
        logical,  optional, intent(in) :: cal

        character(len=56) :: tstep_method
        real(wp) :: tstep_const, time_init, time_end, tref
        logical  :: is_cal

        tref = 1950.0_wp
        if (present(time_ref)) tref = time_ref
        is_cal = .false.
        if (present(cal)) is_cal = cal

        call nml_read(path_par, group, "tstep_method", tstep_method)
        call nml_read(path_par, group, "tstep_const",  tstep_const)
        call nml_read(path_par, group, "time_init",    time_init)
        call nml_read(path_par, group, "time_end",     time_end)
        call nml_read(path_par, group, "dtt",          dtt)

        if (is_cal) then
            call tstep_init(ts, time_init, time_end, method=tstep_method, units="year", &
                            time_ref=tref, const_rel=0.0_wp, const_cal=tstep_const)
        else
            call tstep_init(ts, time_init, time_end, method=tstep_method, units="year", &
                            time_ref=tref, const_rel=tstep_const)
        end if
    end subroutine timeline_init

    function tstep_due(time, dt) result(due)
        ! Cadence predicate: true when `time` falls on the dt grid (0.01-yr
        ! precision). dt <= 0 disables the cadence (never due).
        real(wp), intent(in) :: time, dt
        logical :: due
        due = .false.
        if (dt > 0.0_wp) due = (mod(nint(time*100), nint(dt*100)) == 0)
    end function tstep_due

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

    subroutine domain_startup(dom, ts, bsl, restore_bsl)
        ! Establish the domain state after domain_init: cold start (ctl%restart
        ! == "None") builds the initial boundary state; otherwise the restart
        ! bundle is restored and the hi-res hub rebuilt from the restored models.
        ! restore_bsl (default .true.) also restores the shared bsl from the same
        ! bundle folder -- the single-domain convention, where the run-level
        ! bsl_restart.nc lives in the domain's bundle. Multi-domain drivers
        ! restore the bsl once themselves (bsl_startup) and pass .false..
        ! Flavor drivers with their own cold start (esm, rembo) keep their own
        ! cold branch and call this for the restart branch only.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        type(bsl_class),   intent(inout) :: bsl
        logical, intent(in), optional    :: restore_bsl

        logical :: do_bsl

        do_bsl = .true.
        if (present(restore_bsl)) do_bsl = restore_bsl

        if (trim(dom%ctl%restart) == "None") then
            call domain_init_state(dom, ts, bsl)
        else
            if (do_bsl) call bsl_startup(bsl, ts, trim(dom%ctl%restart))
            call domain_restart_read(dom, trim(dom%ctl%restart), ts, bsl)
            call refresh_htopo(dom)
        end if
    end subroutine domain_startup

    subroutine run_restart_write(dom, bsl, time)
        ! Single-domain restart: write the domain bundle and the run-level shared
        ! bsl restart into the same auto-named per-time folder. Multi-domain
        ! drivers write per-domain bundles + one run-root bsl bundle themselves.
        type(ice_domain), intent(inout) :: dom
        type(bsl_class),  intent(inout) :: bsl
        real(wp),         intent(in)    :: time

        call domain_restart_write(dom, time)
        call bsl_restart_write(bsl, trim(restart_bundle_dir(time))//"/bsl_restart.nc", time)
    end subroutine run_restart_write

    subroutine domain_init(dom, path_par, time, group_suffix, init_climate, timeline_group)
        ! Initialize all sub-models of one domain, load the hi-res reference hub,
        ! prime the Yelmo<->hub maps, and place marine_shelf on its configured grid.
        ! The barystatic sea level (bsl) is NOT a domain sub-model: it is a shared,
        ! driver-owned object (one per run, common to every domain), so it is
        ! initialized by the driver and passed into the isostasy steps.
        !
        ! group_suffix (optional, default "") is appended to every namelist group
        ! name (yelmo -> yelmo<suffix>, coupling -> coupling<suffix>, ...), so
        ! several domains can share one parameter file with disjoint group names
        ! (the multi-domain / bipolar convention). Yelmo physics sub-groups (ydyn,
        ! ytopo, ...) stay shared: they are named by pointer fields inside the
        ! [yelmo<suffix>] block, so the nml decides whether they are shared.
        !
        ! init_climate (optional, default .TRUE.) initializes the snapclim climate
        ! sub-model. Variants that supply their own climate forcing (e.g. the ESM
        ! driver, which owns an esm_forcing_class in place of dom%snp) pass .FALSE.
        ! to skip snapclim_init; grid_clim is still resolved so grid_smb can default
        ! to it.
        !
        ! timeline_group (optional, default "ctrl") names the group holding the
        ! run's shared timeline -- the same group the driver passes to
        ! timeline_init -- from which the domain reads tstep_method/dtt itself.
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: path_par
        real(wp),         intent(in)    :: time       ! model time
        character(len=*), intent(in), optional :: group_suffix
        logical,          intent(in), optional :: init_climate
        character(len=*), intent(in), optional :: timeline_group

        character(len=256)    :: domain, tgroup
        character(len=64)     :: sfx
        logical               :: do_climate
        type(grid_class)      :: grid_m, grid_y, grid_i, grid_c, grid_s
        integer               :: nx_m, ny_m, nx_i, ny_i, nx_c, ny_c, nx_s, ny_s
        real(wp), allocatable :: regions_m(:,:), basins_m(:,:), basins_c(:,:)
        real(wp), allocatable :: xs(:), ys(:), lats_s(:,:), Href_s(:,:)

        sfx = ""
        if (present(group_suffix)) sfx = trim(group_suffix)

        do_climate = .TRUE.
        if (present(init_climate)) do_climate = init_climate

        tgroup = "ctrl"
        if (present(timeline_group)) tgroup = trim(timeline_group)

        ! --- run control ---
        call domain_ctl_load(dom%ctl, path_par, trim(sfx), trim(tgroup))

        ! --- ice sheet (grid read from file) ---
        call yelmo_init(dom%yelmo, filename=path_par, grid_def="file", time=time, &
                        group="yelmo"//trim(sfx))
        domain = trim(dom%yelmo%par%domain)
        dom%ctl%domain     = trim(domain)
        dom%ctl%grid_yelmo = trim(dom%yelmo%par%grid_name)

        ! --- external forcing models (climate/smb/isostasy on the Yelmo grid) ---
        ! Isostasy on its configured grid ([coupling] grid_isos; default = grid_yelmo).
        if (len_trim(dom%ctl%grid_isos) == 0) dom%ctl%grid_isos = trim(dom%ctl%grid_yelmo)
        call grid_cdo_read_desc(grid_i, trim(dom%ctl%grid_isos),  MAP_FLDR)
        call grid_cdo_read_desc(grid_y, trim(dom%ctl%grid_yelmo), MAP_FLDR)
        nx_i = grid_i%G%nx
        ny_i = grid_i%G%ny
        ! Grid spacing in Yelmo units, scaled by the resolution ratio (per axis).
        dom%ctl%dx_isos = dom%yelmo%grd%G%dx * (grid_i%G%dx / grid_y%G%dx)
        dom%ctl%dy_isos = dom%yelmo%grd%G%dy * (grid_i%G%dy / grid_y%G%dy)
        call isos_init(dom%isos, path_par, "isos"//trim(sfx), nx_i, ny_i, &
                       dom%ctl%dx_isos, dom%ctl%dy_isos)

        call sediments_init(dom%sed, path_par, dom%yelmo%grd%G%nx, dom%yelmo%grd%G%ny, &
                            domain, dom%yelmo%par%grid_name, group="sed"//trim(sfx))
        dom%yelmo%bnd%H_sed = dom%sed%now%H

        call geothermal_init(dom%gthrm, path_par, dom%yelmo%grd%G%nx, dom%yelmo%grd%G%ny, &
                             domain, dom%yelmo%par%grid_name, group="ghf"//trim(sfx))
        dom%yelmo%bnd%Q_geo = dom%gthrm%now%ghf

        ! --- hi-res reference hub + coupler ---
        call htopo_init(dom%topo, path_par, "htopo"//trim(sfx))
        dom%ctl%grid_name = trim(dom%topo%par%grid_name)
        if (len_trim(dom%ctl%grid_mshlf) == 0) dom%ctl%grid_mshlf = trim(dom%ctl%grid_name)

        ! Grids resolve from maps/grid_<name>.txt; prime the Yelmo<->hub maps.
        call coupler_init(dom%cpl)
        call coupler_prime(dom%cpl, dom%ctl%grid_yelmo, dom%ctl%grid_name, "bilin")  ! Yelmo -> hub
        call coupler_prime(dom%cpl, dom%ctl%grid_name, dom%ctl%grid_yelmo, "con")    ! hub -> Yelmo

        ! --- climate on its configured grid ([coupling] grid_clim; default = grid_yelmo) ---
        ! snapclim reads grid-specific input data, so grid_clim must be a grid whose
        ! forcing files exist (the Yelmo grid for the standard setup).
        if (len_trim(dom%ctl%grid_clim) == 0) dom%ctl%grid_clim = trim(dom%ctl%grid_yelmo)
        call grid_cdo_read_desc(grid_c, trim(dom%ctl%grid_clim), MAP_FLDR)
        nx_c = grid_c%G%nx
        ny_c = grid_c%G%ny
        dom%ctl%dx_clim = dom%yelmo%grd%G%dx * (grid_c%G%dx / grid_y%G%dx)
        if (do_climate) then
            call remap(dom, dom%topo%basins, dom%ctl%grid_name, basins_c, dom%ctl%grid_clim, "nn")
            call snapclim_init(dom%snp, path_par, domain, trim(dom%ctl%grid_clim), &
                               nx_c, ny_c, basins_c, group="snap"//trim(sfx))
        end if

        ! --- smb on its configured grid ([coupling] grid_smb; default = grid_clim) ---
        ! smbpal reads no grid-specific data; only lats (insolation) is physical.
        if (len_trim(dom%ctl%grid_smb) == 0) dom%ctl%grid_smb = trim(dom%ctl%grid_clim)
        call grid_cdo_read_desc(grid_s, trim(dom%ctl%grid_smb), MAP_FLDR)
        nx_s = grid_s%G%nx
        ny_s = grid_s%G%ny
        allocate(xs(nx_s), ys(ny_s), lats_s(nx_s, ny_s))
        xs     = real(grid_s%G%x, wp)
        ys     = real(grid_s%G%y, wp)
        lats_s = real(grid_s%lat, wp)
        call smbpal_init(dom%smb, path_par, x=xs, y=ys, lats=lats_s, &
                         group="smbpal"//trim(sfx), itm_group="itm"//trim(sfx))

        ! Alternative SMB (smb_simple) on the same grid, if selected. Unlike
        ! smbpal (1D axes), smb_simple takes 2D projected coordinates.
        if (trim(dom%ctl%smb_method) == "smb_simple") then
            call smb_simple_init(dom%smbs, path_par, x=real(grid_s%x, wp), &
                                 y=real(grid_s%y, wp), lat=lats_s, &
                                 group="smb_simple"//trim(sfx), units="m")
            call remap(dom, dom%yelmo%bnd%H_ice_ref, dom%ctl%grid_yelmo, &
                       Href_s, dom%ctl%grid_smb, "bilin")
            call smb_simple_set_mask(dom%smbs, Href_s)
        end if

        ! --- marine_shelf on its configured grid (grid_y already read above) ---
        call grid_cdo_read_desc(grid_m, trim(dom%ctl%grid_mshlf), MAP_FLDR)
        nx_m = grid_m%G%nx
        ny_m = grid_m%G%ny
        ! Grid spacing in Yelmo dx units, scaled by the resolution ratio.
        dom%ctl%dx_mshlf = dom%yelmo%grd%G%dx * (grid_m%G%dx / grid_y%G%dx)

        ! Region/basin masks on the mshlf grid (from the hub).
        call remap(dom, dom%topo%regions, dom%ctl%grid_name, regions_m, dom%ctl%grid_mshlf, "nn")
        call remap(dom, dom%topo%basins,  dom%ctl%grid_name, basins_m,  dom%ctl%grid_mshlf, "nn")

        call marshelf_init(dom%mshlf, path_par, "marine_shelf"//trim(sfx), nx_m, ny_m, &
                           domain, trim(dom%ctl%grid_mshlf), regions_m, basins_m)

        ! Optimization state (basal friction + thermal forcing); no-op unless
        ! equil_method == "opt". Must follow yelmo_init (grid + till params known).
        call domain_opt_init(dom, path_par, trim(sfx))

        ! NEGIS cb_ref modification (Greenland): load its [negis] parameters when
        ! enabled, so use_negis=True cannot silently run with default factors.
        if (dom%ctl%use_negis) call negis_par_load(dom%ngs, path_par, trim(sfx))

    end subroutine domain_init

    subroutine negis_par_load(ngs, path_par, suffix)
        ! Load the NEGIS cb_ref scaling parameters ([negis<suffix>]). Ported from
        ! yelmox.f90; only read when [coupling] use_negis is set.
        type(negis_params), intent(inout) :: ngs
        character(len=*),   intent(in)    :: path_par
        character(len=*),   intent(in)    :: suffix

        ngs%use_negis_par = .true.
        call nml_read(path_par, "negis"//trim(suffix), "cf_0",      ngs%cf_0)
        call nml_read(path_par, "negis"//trim(suffix), "cf_1",      ngs%cf_1)
        call nml_read(path_par, "negis"//trim(suffix), "cf_centre", ngs%cf_centre)
        call nml_read(path_par, "negis"//trim(suffix), "cf_north",  ngs%cf_north)
        call nml_read(path_par, "negis"//trim(suffix), "cf_south",  ngs%cf_south)
    end subroutine negis_par_load

    subroutine domain_opt_init(dom, path_par, suffix)
        ! Load optimization parameters and prepare Yelmo for external cb_ref:
        ! allocate/seed the friction bounds (cf_min/cf_max) on the Yelmo grid and
        ! switch till_method to external (-1) so yelmo_update uses the optimized
        ! cb_ref. The initial cb_ref guess (cold start) is set in domain_init_state
        ! after yelmo_init_state; on restart cb_ref is restored from the bundle.
        ! No-op unless equil_method == "opt".
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: path_par
        character(len=*), intent(in)    :: suffix

        integer :: nx, ny

        if (trim(dom%ctl%equil_method) /= "opt") return

        dom%opt%tf_basins = 0
        call optimize_par_load(dom%opt, path_par, "opt"//trim(suffix))

        nx = dom%yelmo%grd%G%nx
        ny = dom%yelmo%grd%G%ny
        allocate(dom%opt%cf_min(nx, ny), dom%opt%cf_max(nx, ny))
        dom%opt%cf_min = dom%yelmo%dyn%par%till_cf_min
        dom%opt%cf_max = dom%yelmo%dyn%par%till_cf_ref

        ! cb_ref is set externally by the optimization from here on.
        dom%yelmo%dyn%par%till_method = -1
    end subroutine domain_opt_init

    subroutine domain_regions_init(dom, outfldr)
        ! Define the domain's regions of interest for 1D regional output. Masks are
        ! resolved on the Yelmo grid (get_ice_sub_region); regional files land in
        ! outfldr. Domains without defined sub-regions get n=0 (global region only).
        ! Must be called after domain_init and before the first yelmo_update.
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr

        logical, allocatable :: tmp_mask(:,:)
        character(len=256)   :: domain, grid_name
        integer              :: i, nx, ny

        domain    = trim(dom%ctl%domain)
        grid_name = trim(dom%ctl%grid_yelmo)
        nx = dom%yelmo%grd%G%nx
        ny = dom%yelmo%grd%G%ny
        allocate(tmp_mask(nx, ny))

        select case(trim(domain))

            case("Antarctica")
                call yelmo_regions_init(dom%yelmo, n=3)

                call get_ice_sub_region(tmp_mask, "APIS", domain, grid_name)
                call yelmo_region_init(dom%yelmo%regs(1), "APIS", mask=tmp_mask, &
                                       write_to_file=.true., outfldr=outfldr)

                call get_ice_sub_region(tmp_mask, "WAIS", domain, grid_name)
                call yelmo_region_init(dom%yelmo%regs(2), "WAIS", mask=tmp_mask, &
                                       write_to_file=.true., outfldr=outfldr)

                call get_ice_sub_region(tmp_mask, "EAIS", domain, grid_name)
                call yelmo_region_init(dom%yelmo%regs(3), "EAIS", mask=tmp_mask, &
                                       write_to_file=.true., outfldr=outfldr)

            case("Laurentide")
                ! Prevent ice growth in Greenland (region 1.30) and on grid borders.
                where(abs(dom%yelmo%bnd%regions - 1.30) < 1e-3) &
                    dom%yelmo%bnd%mask_ice = MASK_ICE_NONE
                dom%yelmo%bnd%mask_ice(1,:)  = MASK_ICE_NONE
                dom%yelmo%bnd%mask_ice(nx,:) = MASK_ICE_NONE
                dom%yelmo%bnd%mask_ice(:,1)  = MASK_ICE_NONE
                dom%yelmo%bnd%mask_ice(:,ny) = MASK_ICE_NONE

                call yelmo_regions_init(dom%yelmo, n=1)
                call get_ice_sub_region(tmp_mask, "Hudson", domain, grid_name)
                call yelmo_region_init(dom%yelmo%regs(1), "Hudson", mask=tmp_mask, &
                                       write_to_file=.true., outfldr=outfldr)

            case("Greenland")
                ! Prevent ice in Iceland/Svalbard (regions 1.20/1.23/1.31, grid borders).
                where(abs(dom%yelmo%bnd%regions - 1.20) < 1e-3) dom%yelmo%bnd%mask_ice = MASK_ICE_NONE
                where(abs(dom%yelmo%bnd%regions - 1.23) < 1e-3) dom%yelmo%bnd%mask_ice = MASK_ICE_NONE
                where(abs(dom%yelmo%bnd%regions - 1.31) < 1e-3) dom%yelmo%bnd%mask_ice = MASK_ICE_NONE

                ! NEGIS cb_ref modification: enabled via [coupling] use_negis, which
                ! loads the [negis] parameters in domain_init.

                ! With external cb_ref (till_method=-1) start from the reference value.
                if (dom%yelmo%dyn%par%till_method == -1) &
                    dom%yelmo%dyn%now%cb_ref = dom%yelmo%dyn%par%till_cf_ref

                call yelmo_regions_init(dom%yelmo, n=0)

            case("Patagonia")
                ! Fix ice on the domain borders, relax to obs outside the icefield.
                dom%yelmo%bnd%mask_ice        = MASK_ICE_DYNAMIC
                dom%yelmo%bnd%mask_ice(1,:)   = MASK_ICE_FIXED
                dom%yelmo%bnd%mask_ice(nx,:)  = MASK_ICE_FIXED
                dom%yelmo%bnd%mask_ice(:,1)   = MASK_ICE_FIXED
                dom%yelmo%bnd%mask_ice(:,ny)  = MASK_ICE_FIXED
                where(abs(dom%yelmo%bnd%regions - 1.0) < 1e-3)
                    dom%yelmo%bnd%tau_relax = -1.0      ! icefield: free evolution
                elsewhere
                    dom%yelmo%bnd%tau_relax = 50.0      ! outside: relax to H_ice_ref
                end where
                call yelmo_regions_init(dom%yelmo, n=0)

            case default
                ! No sub-regions defined for this domain; global region only.
                call yelmo_regions_init(dom%yelmo, n=0)

        end select

        ! Rename the regional 1D files to the climber-x timeseries convention:
        !   global -> yelmo_<grid>_ts.nc, sub-region k -> yelmo_<grid>_ts_<name>.nc
        dom%yelmo%reg%fnm = trim(outfldr)//"yelmo_"//trim(grid_name)//"_ts.nc"
        if (dom%yelmo%par%n_reg > 0) then
            do i = 1, dom%yelmo%par%n_reg
                dom%yelmo%regs(i)%fnm = trim(outfldr)//"yelmo_"//trim(grid_name)// &
                                        "_ts_"//trim(dom%yelmo%regs(i)%name)//".nc"
            end do
        end if

    end subroutine domain_regions_init

    subroutine domain_init_state(dom, ts, bsl)
        ! Build the initial boundary state and initialize the Yelmo state
        ! variables. Lifted from yelmox.f90's "update initial boundary conditions"
        ! block (minimal core: no smb_simple / domain-special / optimization).
        ! bsl is the shared, driver-owned sea level; the driver has already called
        ! bsl_update for the initial time, so this routine only consumes it.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        type(bsl_class),   intent(inout) :: bsl

        real(wp), allocatable :: z_bed_ref_i(:,:), H_ice_ref_i(:,:)
        real(wp), allocatable :: z_bed_i(:,:), H_ice_i(:,:)
        real(wp), allocatable :: z_srf_c(:,:), basins_c(:,:)
        character(len=256) :: gi, gy, gc, gn

        gi = trim(dom%ctl%grid_isos)
        gy = trim(dom%ctl%grid_yelmo)
        gc = trim(dom%ctl%grid_clim)
        gn = trim(dom%ctl%grid_name)

        ! Sea level + isostasy reference state (isostasy runs on grid_isos)
        call remap(dom, dom%yelmo%bnd%z_bed_ref, gy, z_bed_ref_i, gi, "bilin")
        call remap(dom, dom%yelmo%bnd%H_ice_ref, gy, H_ice_ref_i, gi, "bilin")
        call isos_init_ref(dom%isos, z_bed_ref_i, H_ice_ref_i)
        call remap(dom, dom%yelmo%bnd%z_bed,      gy, z_bed_i,     gi, "bilin")
        call remap(dom, dom%yelmo%tpo%now%H_ice,  gy, H_ice_i,     gi, "bilin")
        call isos_init_state(dom%isos, z_bed_i, H_ice_i, ts%time, bsl)
        call couple_isostasy_to_yelmo(dom)

        ! Refresh the hub from the initial geometry; climate/smb/mshlf read from it.
        call refresh_htopo(dom)

        ! Climate on grid_clim (note: init uses time_rel for snapclim), then the
        ! surface mass balance on grid_smb (smbpal or smb_simple; init=.true.
        ! runs the smbpal ITM equilibration before the first update).
        if (dom%ctl%with_climate) then
            call remap(dom, dom%topo%z_srf,  gn, z_srf_c,  gc, "bilin")
            call remap(dom, dom%topo%basins, gn, basins_c, gc, "nn")
            call snapclim_update(dom%snp, z_srf=z_srf_c, time=ts%time_rel, &
                                 domain=dom%ctl%domain, dx=dom%ctl%dx_clim, basins=basins_c)
            call domain_update_smb(dom, ts, init=.true.)
        end if

        ! Marine shelf through the (already refreshed) hub.
        call step_marine_shelf(dom, ts)

        ! Assemble the Yelmo boundary state from the freshly produced module
        ! outputs (smb + marine_shelf; isostasy already coupled above).
        call couple_smb_to_yelmo(dom)
        call couple_marine_to_yelmo(dom)

        ! Initialize state variables (dyn, therm, mat) with a cold base
        call yelmo_print_bound(dom%yelmo%bnd)
        call yelmo_init_state(dom%yelmo, time=ts%time, thrm_method="robin-cold")

        ! Cold-start friction guess for the optimization (restart restores cb_ref).
        if (trim(dom%ctl%equil_method) == "opt") then
            dom%yelmo%dyn%now%cb_ref = dom%opt%cf_init
        end if

        ! Domain-specific cold-start setup (equilibration / LGM initialization /
        ! Greenland marine-ice). Cold start only; restart skips it.
        call domain_init_special(dom, ts)

    end subroutine domain_init_state

    subroutine domain_init_special(dom, ts)
        ! Domain-specific cold-start startup, dispatched on domain name. Mirrors
        ! yelmox.f90's per-domain startup block. The DEFAULT (incl. Antarctica)
        ! path runs a short equilibration to synchronize the model fields.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        select case(trim(dom%ctl%domain))

            case("Laurentide")
                ! Steady-state: LGM reconstruction; transient: grow from zero ice.
                if (trim(dom%ctl%tstep_method) == "const") then
                    call domain_init_lgm_north(dom, ts, "Laurentide", "ref_lgm")
                else
                    call domain_init_lgm_north(dom, ts, "Laurentide", "zero")
                end if

            case("North")
                ! Steady-state only: whole-NH LGM reconstruction (ICE-6G_C).
                if (trim(dom%ctl%tstep_method) == "const") then
                    call domain_init_lgm_north(dom, ts, "North", "ref_lgm")
                end if

            case("Greenland")
                ! Optionally impose LGM-like marine ice; otherwise no startup equil.
                if (dom%ctl%greenland_init_marine_H) then
                    where(dom%yelmo%bnd%mask_ice /= MASK_ICE_NONE .and. &
                          dom%yelmo%tpo%now%H_ice < 600.0_wp .and. &
                          dom%yelmo%bnd%z_bed > -500.0_wp)
                        dom%yelmo%tpo%now%H_ice = 800.0_wp
                    end where
                    if (dom%ctl%with_ice_sheet) &
                        call yelmo_update_equil(dom%yelmo, ts%time, time_tot=10.0_wp, &
                                                dt=1.0_wp, topo_fixed=.FALSE.)
                end if

            case default
                ! Antarctica etc.: short equilibration with constant boundaries.
                if (dom%ctl%with_ice_sheet) &
                    call yelmo_update_equil(dom%yelmo, ts%time, time_tot=10.0_wp, &
                                            dt=1.0_wp, topo_fixed=.FALSE.)

        end select

    end subroutine domain_init_special

    subroutine domain_init_lgm_north(dom, ts, region, method)
        ! Initialize a Northern-Hemisphere domain (Laurentide or whole "North")
        ! from the ICE-6G_C LGM reconstruction. Sets the reconstructed grounded
        ! ice as the initial thickness (method-dependent), refreshes the surface
        ! and (via the hub) the climate/smb, and stabilizes the dynamic fields.
        ! Ported from yelmox.f90's yelmox_init_{laurentide,north}_lgm.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        character(len=*),  intent(in)    :: region   ! "Laurentide" or "North"
        character(len=*),  intent(in)    :: method   ! "const", "ref_lgm", else zero

        character(len=1024) :: path_lgm, grid_name
        integer  :: nx, ny
        real(wp) :: beta_min_save

        nx = dom%yelmo%tpo%par%nx
        ny = dom%yelmo%tpo%par%ny
        grid_name = trim(dom%yelmo%par%grid_name)

        ! Load LGM reconstruction (slice 1) into the reference ice thickness.
        path_lgm = "ice_data/"//trim(region)//"/"//trim(grid_name)//"/"// &
                   trim(grid_name)//"_TOPO-ICE-6G_C.nc"
        call nc_read(path_lgm, "dz", dom%yelmo%bnd%H_ice_ref, start=[1,1,1], &
                     count=[nx,ny,1])

        ! Determine the initial ice thickness.
        select case(trim(method))
            case("const")
                dom%yelmo%tpo%now%H_ice = 0.0_wp
                where (dom%yelmo%bnd%regions == 1.1_wp .and. dom%yelmo%bnd%z_bed > 0.0_wp) &
                    dom%yelmo%tpo%now%H_ice = 1000.0_wp
                where (dom%yelmo%bnd%regions == 1.12_wp) dom%yelmo%tpo%now%H_ice = 1000.0_wp
                call smooth_gauss_2D(dom%yelmo%tpo%now%H_ice, dx=real(dom%yelmo%grd%G%dx,wp), f_sigma=3.0_wp)
                call yelmo_init_topo(dom%yelmo, trim(dom%ctl%path_par), &
                                     dom%yelmo%par%nml_init_topo, ts%time, load_topo=.FALSE.)
            case("ref_lgm")
                where ( dom%yelmo%bnd%z_bed > -500.0_wp .and. &
                        (dom%yelmo%bnd%regions == 1.1_wp  .or. &
                         dom%yelmo%bnd%regions == 1.11_wp .or. &
                         dom%yelmo%bnd%regions == 1.12_wp) )
                    dom%yelmo%tpo%now%H_ice = dom%yelmo%bnd%H_ice_ref
                end where
                call smooth_gauss_2D(dom%yelmo%tpo%now%H_ice, dx=real(dom%yelmo%grd%G%dx,wp), f_sigma=2.0_wp)
                call yelmo_init_topo(dom%yelmo, trim(dom%ctl%path_par), &
                                     dom%yelmo%par%nml_init_topo, ts%time, load_topo=.FALSE.)
            case default
                ! Zero ice thickness (transient start): do nothing.
        end select

        ! Update surface topography fields (fixed H), then remove thin floating ice.
        call yelmo_update_equil(dom%yelmo, ts%time, time_tot=1.0_wp, dt=1.0_wp, topo_fixed=.TRUE.)
        where(dom%yelmo%tpo%now%mask_bed == 5 .and. dom%yelmo%tpo%now%H_ice < 50.0_wp) &
            dom%yelmo%tpo%now%H_ice = 0.0_wp
        call yelmo_update_equil(dom%yelmo, ts%time, time_tot=1.0_wp, dt=1.0_wp, topo_fixed=.TRUE.)

        if (trim(method) == "ref_lgm") then
            ! Store the clean thickness as the reference state (drives smb masks).
            dom%yelmo%bnd%H_ice_ref = dom%yelmo%tpo%now%H_ice
        end if

        ! Refresh the hub and climate/smb to reflect the new geometry, then land
        ! the smb on the Yelmo grid (used/adjusted just below).
        call refresh_htopo(dom)
        call step_climate(dom, ts)
        call couple_smb_to_yelmo(dom)

        if (trim(method) == "const") then
            ! Ensure ice can grow on high-latitude land (mainly Cordilleran).
            where (dom%yelmo%bnd%regions == 1.1_wp .and. dom%yelmo%grd%lat > 50.0_wp .and. &
                   dom%yelmo%bnd%z_bed > 0.0_wp .and. dom%yelmo%bnd%smb < 0.0_wp) &
                dom%yelmo%bnd%smb = 0.5_wp
            if (dom%ctl%with_ice_sheet) &
                call yelmo_update_equil(dom%yelmo, ts%time, time_tot=5e3_wp, dt=5.0_wp, &
                                        topo_fixed=.FALSE.)
        else
            ! ref_lgm / zero: stabilize dynamic fields with a raised beta_min.
            if (dom%ctl%with_ice_sheet) then
                beta_min_save = dom%yelmo%dyn%par%beta_min
                dom%yelmo%dyn%par%beta_min = 100.0_wp
                call yelmo_update_equil(dom%yelmo, ts%time, time_tot=2e2_wp, dt=5.0_wp, &
                                        topo_fixed=.FALSE.)
                dom%yelmo%dyn%par%beta_min = beta_min_save
            end if
        end if

    end subroutine domain_init_lgm_north

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
        ! not written -- it is rebuilt by refresh_htopo from the restored models.
        ! The shared barystatic sea level is NOT written here -- the driver owns it
        ! and writes a single bsl_restart.nc for the whole run.
        ! `outfldr` (optional) prefixes the auto-named per-time folder, so each
        ! domain of a multi-domain run writes into its own subfolder.
        type(ice_domain), intent(inout) :: dom
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
        type(ice_domain),  intent(inout) :: dom
        character(len=*),  intent(in)    :: fldr
        type(tstep_class), intent(in)    :: ts
        type(bsl_class),   intent(inout) :: bsl

        real(wp), allocatable :: z_bed_i(:,:), H_ice_i(:,:)
        character(len=256) :: gi, gy

        gi = trim(dom%ctl%grid_isos)
        gy = trim(dom%ctl%grid_yelmo)

        ! Restore Yelmo first: it provides the current H_ice/z_bed for isostasy.
        ! Two reads are needed, mirroring yelmo's native init-from-restart:
        !   - yelmo_restart_read_topo_bnd loads the geometry [tpo]+[bnd]
        !     (H_ice, z_bed, ...); the standalone yelmo_restart_read does NOT.
        !   - yelmo_restart_read loads [dyn,therm,mat] + mask_bed.
        ! use_restart/pc_active are flags the native path sets; the topo
        ! diagnostics (f_ice/f_grnd/H_grnd/z_srf) are reconciled below.
        call yelmo_restart_read_topo_bnd(dom%yelmo%tpo, dom%yelmo%bnd, dom%yelmo%time, &
                dom%yelmo%par%restart_interpolated, dom%yelmo%grd, dom%yelmo%par%domain, &
                dom%yelmo%par%grid_name, trim(fldr)//"/yelmo_restart.nc", ts%time)
        call yelmo_restart_read(dom%yelmo, trim(fldr)//"/yelmo_restart.nc", ts%time)
        dom%yelmo%par%use_restart = .true.
        dom%yelmo%time%pc_active  = .true.

        ! Restore isostasy via isos_init_state (reads state + reference from the
        ! bundle and runs the full post-read setup), on the isos grid. The shared
        ! bsl was already restored + updated by the driver before this call.
        dom%isos%par%use_restart = .true.
        dom%isos%par%restart     = trim(fldr)//"/isos_restart.nc"
        call remap(dom, dom%yelmo%bnd%z_bed,     gy, z_bed_i, gi, "bilin")
        call remap(dom, dom%yelmo%tpo%now%H_ice, gy, H_ice_i, gi, "bilin")
        call isos_init_state(dom%isos, z_bed_i, H_ice_i, ts%time, bsl)
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

    subroutine domain_ctl_load(ctl, path_par, suffix, timeline_group)
        ! Load this domain's setup + coupling + output config. All groups carry an
        ! optional domain suffix (e.g. "_north"), so several domains can coexist in
        ! one parameter file without group.name collisions (matters for runme -p).
        ! The shared timeline is driver-owned (timeline_init); the values the
        ! domain logic needs (tstep_method, dtt) are read from the same
        ! timeline_group here, so nothing is injected after init.
        type(domain_ctl), intent(inout) :: ctl
        character(len=*), intent(in)    :: path_par
        character(len=*), intent(in)    :: suffix
        character(len=*), intent(in)    :: timeline_group

        character(len=256) :: gc, go

        ctl%path_par = trim(path_par)

        ! Shared timeline values used by the domain logic (dt_clim cadence,
        ! optimization dt, domain-specific startup).
        call nml_read(path_par, timeline_group, "tstep_method", ctl%tstep_method)
        call nml_read(path_par, timeline_group, "dtt",          ctl%dtt)

        ! Domain setup + coupling ([coupling<suffix>]): active components, methods,
        ! per-component grids, restart bundle, restart cadence.
        gc = "coupling"//trim(suffix)
        call nml_read(path_par, gc, "with_ice_sheet",    ctl%with_ice_sheet)
        call nml_read(path_par, gc, "with_isostasy",     ctl%with_isostasy)
        call nml_read(path_par, gc, "with_climate",      ctl%with_climate)
        call nml_read(path_par, gc, "with_marine_shelf", ctl%with_marine_shelf)
        call nml_read(path_par, gc, "equil_method",   ctl%equil_method)
        ctl%smb_method = "smbpal"
        call nml_read(path_par, gc, "smb_method",     ctl%smb_method)
        call nml_read(path_par, gc, "dt_restart",     ctl%dt_restart)
        call nml_read(path_par, gc, "dt_clim",        ctl%dt_clim)

        ! Domain-specific startup / physics switches (Greenland only; keep False
        ! elsewhere). use_negis additionally requires a [negis<suffix>] group.
        call nml_read(path_par, gc, "scale_glacial_smb",       ctl%scale_glacial_smb)
        call nml_read(path_par, gc, "lim_pd_ice",              ctl%lim_pd_ice)
        call nml_read(path_par, gc, "use_negis",               ctl%use_negis)
        call nml_read(path_par, gc, "greenland_init_marine_H", ctl%greenland_init_marine_H)

        ctl%grid_mshlf = ""
        call nml_read(path_par, gc, "grid_mshlf",     ctl%grid_mshlf)
        ctl%grid_isos = ""
        call nml_read(path_par, gc, "grid_isos",      ctl%grid_isos)
        ctl%grid_clim = ""
        call nml_read(path_par, gc, "grid_clim",      ctl%grid_clim)
        ctl%grid_smb = ""
        call nml_read(path_par, gc, "grid_smb",       ctl%grid_smb)
        ctl%restart = "None"
        call nml_read(path_par, gc, "restart",        ctl%restart)

        ! Per-module output switches ([output<suffix>]); default = write everything.
        go = "output"//trim(suffix)
        call nml_read(path_par, go, "write_yelmo", ctl%write_yelmo)
        call nml_read(path_par, go, "write_isos",  ctl%write_isos)
        call nml_read(path_par, go, "write_mshlf", ctl%write_mshlf)
        call nml_read(path_par, go, "write_smb",   ctl%write_smb)
        call nml_read(path_par, go, "write_snap",  ctl%write_snap)
        call nml_read(path_par, go, "write_htopo", ctl%write_htopo)
    end subroutine domain_ctl_load

    subroutine yelmox_step(dom, ts, bsl)
        ! Advance one domain by one coupling step. Component order is fixed here;
        ! each step_* is a no-op when its ctl flag is off. bsl is the shared,
        ! driver-owned sea level (the driver calls bsl_update once per step before
        ! this call); step_isostasy consumes it. Used by the single-domain driver;
        ! the multi-domain driver interleaves the step_* primitives itself.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        type(bsl_class),   intent(inout) :: bsl

        call step_optimize(dom, ts)      ! spinup relaxation + cb_ref/tf_corr tuning
        call step_isostasy(dom, ts, bsl)
        call step_icesheet(dom, ts)
        call refresh_htopo(dom)          ! hi-res geometry mirror, from the models
        call step_climate(dom, ts)       ! climate/smb read geometry from the hub
        call step_marine_shelf(dom, ts)
    end subroutine yelmox_step

    subroutine step_optimize(dom, ts)
        ! Spin-up tuning (equil_method == "opt"): ramp the topography relaxation
        ! timescale, then nudge the basal-friction field cb_ref and the marine
        ! thermal-forcing correction tf_corr toward present-day observations.
        !
        ! cb_ref is a Yelmo-grid control, optimized in place. tf_corr lives on the
        ! marine_shelf grid; the observational targets (H_ice/H_grnd) live on the
        ! Yelmo grid, so the correction is lifted to the Yelmo grid (tf_corr_y),
        ! optimized there, and remapped back to the shelf grid. At identity grids
        ! both remaps are copies, reproducing yelmox.f90 exactly.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        real(wp), allocatable :: tf_corr_y(:,:), tf_corr_m(:,:)
        character(len=256) :: gm, gy

        if (trim(dom%ctl%equil_method) /= "opt") return

        gm = trim(dom%ctl%grid_mshlf)
        gy = trim(dom%ctl%grid_yelmo)

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
    end subroutine step_optimize

    subroutine step_isostasy(dom, ts, bsl)
        ! Run isostasy on its own grid: ice load from Yelmo (bilin). The bedrock /
        ! sea-surface outputs stay on grid_isos (in dom%isos%out); they are landed
        ! on the Yelmo grid by couple_isostasy_to_yelmo (in step_icesheet, before
        ! yelmo_update). Assumes grid_isos is at least as fine as grid_yelmo
        ! (identity when equal). bsl is the shared, driver-owned sea level (already
        ! updated for this step by the driver); isos_update reads it and, under
        ! fastiso/mixed, writes back the prognostic bsl_now -- so with several
        ! domains sharing one bsl the sea level integrates every domain's ice load.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        type(bsl_class),   intent(inout) :: bsl

        real(wp), allocatable :: H_ice_i(:,:), dwdt_i(:,:)
        character(len=256) :: gi, gy

        if (.not. dom%ctl%with_isostasy) return

        gi = trim(dom%ctl%grid_isos)
        gy = trim(dom%ctl%grid_yelmo)

        ! ice load + correction: Yelmo -> isos grid
        call remap(dom, dom%yelmo%tpo%now%H_ice,  gy, H_ice_i, gi, "bilin")
        call remap(dom, dom%yelmo%bnd%dzbdt_corr, gy, dwdt_i,  gi, "bilin")

        call isos_update(dom%isos, H_ice_i, ts%time, bsl, dwdt_corr=dwdt_i)
    end subroutine step_isostasy

    ! --- Yelmo-input couplers -------------------------------------------------
    ! Each coupler remaps one module's output onto the Yelmo grid and assigns it
    ! into yelmo%bnd, i.e. "remap what Yelmo needs, before Yelmo runs". They are
    ! called from step_icesheet (before yelmo_update) and from the init/restart
    ! paths (before yelmo_init_state), so the Yelmo boundary assembly lives in one
    ! place. Each is a no-op when its component is inactive. At identity grids the
    ! remaps are copies, reproducing the pre-refactor behavior exactly.

    subroutine couple_isostasy_to_yelmo(dom)
        ! Bedrock + sea surface from isostasy (grid_isos -> Yelmo, conservative).
        type(ice_domain), intent(inout) :: dom

        real(wp), allocatable :: z_bed_y(:,:), z_ss_y(:,:)
        character(len=256) :: gi, gy

        if (.not. dom%ctl%with_isostasy) return

        gi = trim(dom%ctl%grid_isos)
        gy = trim(dom%ctl%grid_yelmo)

        call remap(dom, dom%isos%out%z_bed, gi, z_bed_y, gy, "con")
        call remap(dom, dom%isos%out%z_ss,  gi, z_ss_y,  gy, "con")
        dom%yelmo%bnd%z_bed = z_bed_y
        dom%yelmo%bnd%z_sl  = z_ss_y
    end subroutine couple_isostasy_to_yelmo

    subroutine couple_smb_to_yelmo(dom)
        ! Surface mass balance + surface temperature from the active SMB model
        ! (grid_smb -> Yelmo, conservative), with the we->ie unit scaling and the
        ! optional Greenland modifications. The producing step (domain_update_smb,
        ! or a flavor climate step) leaves smb/tsrf on grid_smb in the SMB model's
        ! own fields; this coupler is the single place that lands them on Yelmo.
        type(ice_domain), intent(inout) :: dom

        real(wp), allocatable :: smb_y(:,:), tsrf_y(:,:), ta_y(:,:), ta_pd_y(:,:)
        character(len=256) :: gs, gc, gy

        if (.not. dom%ctl%with_climate) return

        gs = trim(dom%ctl%grid_smb)
        gc = trim(dom%ctl%grid_clim)
        gy = trim(dom%ctl%grid_yelmo)

        if (trim(dom%ctl%smb_method) == "smb_simple") then
            call remap(dom, dom%smbs%smb,   gs, smb_y,  gy, "con")
            call remap(dom, dom%smbs%t_srf, gs, tsrf_y, gy, "con")
        else
            call remap(dom, dom%smb%ann%smb,  gs, smb_y,  gy, "con")
            call remap(dom, dom%smb%ann%tsrf, gs, tsrf_y, gy, "con")
        end if

        dom%yelmo%bnd%smb   = smb_y * dom%yelmo%bnd%c%conv_we_ie * 1e-3
        dom%yelmo%bnd%T_srf = tsrf_y

        ! Glacial-smb modification (Greenland): reduce large negative smb toward a
        ! quasi glacial-interglacial index. Operates on the aggregated Yelmo-grid smb.
        if (trim(dom%ctl%domain) == "Greenland" .and. dom%ctl%scale_glacial_smb) then
            call remap(dom, dom%snp%now%ta_ann,   gc, ta_y,    gy, "bilin")
            call remap(dom, dom%snp%clim0%ta_ann, gc, ta_pd_y, gy, "bilin")
            call calc_glacial_smb(dom%yelmo%bnd%smb, real(dom%yelmo%grd%lat,wp), ta_y, ta_pd_y)
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
        type(ice_domain), intent(inout) :: dom

        real(wp), allocatable :: bmb_y(:,:), Tshlf_y(:,:)
        character(len=256) :: gm, gy

        if (.not. dom%ctl%with_marine_shelf) return

        gm = trim(dom%ctl%grid_mshlf)
        gy = trim(dom%ctl%grid_yelmo)

        call remap(dom, dom%mshlf%now%bmb_shlf, gm, bmb_y,   gy, "con")
        call remap(dom, dom%mshlf%now%T_shlf,   gm, Tshlf_y, gy, "con")
        dom%yelmo%bnd%bmb_shlf = bmb_y
        dom%yelmo%bnd%T_shlf   = Tshlf_y
    end subroutine couple_marine_to_yelmo

    subroutine step_icesheet(dom, ts)
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        ! Greenland NEGIS: update cb_ref from bed properties + NEGIS scaling.
        if (trim(dom%ctl%domain) == "Greenland" .and. dom%ngs%use_negis_par) &
            call negis_update_cb_ref(dom%yelmo, dom%ngs, ts%time)

        ! Assemble the Yelmo boundary state: remap every coupled module's output
        ! onto the Yelmo grid, before Yelmo runs. isostasy was produced this step;
        ! smb / marine_shelf were produced last step (the one-step coupling lag).
        call couple_isostasy_to_yelmo(dom)
        call couple_smb_to_yelmo(dom)
        call couple_marine_to_yelmo(dom)

        if (.not. dom%ctl%with_ice_sheet) return
        if (ts%n == 0 .and. dom%yelmo%par%use_restart) return

        call yelmo_update(dom%yelmo, ts%time)
    end subroutine step_icesheet

    subroutine negis_update_cb_ref(ylmo, ngs, time)
        ! Northeast Greenland Ice Stream cb_ref modification: recompute cb_ref from
        ! bed properties (calc_cb_ref), then scale the NEGIS basins (9.1/9.2/9.3)
        ! by time-dependent factors. Ported from yelmox.f90. Requires the [negis]
        ! cf_* parameters to be loaded; disabled by default (see domain_regions_init).
        type(yelmo_class),  intent(inout) :: ylmo
        type(negis_params), intent(inout) :: ngs
        real(wp),           intent(in)    :: time

        integer :: i, j, nx, ny

        nx = ylmo%grd%G%nx
        ny = ylmo%grd%G%ny

        if (time < -11e3_wp) then
            ngs%cf_x = ngs%cf_0
        else
            ngs%cf_x = ngs%cf_0 + (time - (-11e3_wp)) / (0.0_wp - (-11e3_wp)) * (ngs%cf_1 - ngs%cf_0)
        end if

        if (time < -4e3_wp) then
            ngs%cf_south = 1.0_wp
        else
            ngs%cf_north = 1.0_wp
        end if

        ! Recompute cb_ref like the standard till function.
        call calc_cb_ref(ylmo%dyn%now%cb_ref, ylmo%bnd%z_bed, ylmo%bnd%z_bed_sd, ylmo%bnd%z_sl, &
                ylmo%bnd%H_sed, ylmo%dyn%par%till_f_sed, ylmo%dyn%par%till_sed_min, ylmo%dyn%par%till_sed_max, &
                ylmo%dyn%par%till_cf_ref, ylmo%dyn%par%till_cf_min, ylmo%dyn%par%till_z0, ylmo%dyn%par%till_z1, &
                ylmo%dyn%par%till_n_sd, ylmo%dyn%par%till_scale_zb, ylmo%dyn%par%till_scale_sed)

        ! Apply NEGIS basin scaling.
        do j = 1, ny
        do i = 1, nx
            if (ylmo%bnd%basins(i,j) == 9.1_wp) ylmo%dyn%now%cb_ref(i,j) = ylmo%dyn%now%cb_ref(i,j) * ngs%cf_centre
            if (ylmo%bnd%basins(i,j) == 9.2_wp) ylmo%dyn%now%cb_ref(i,j) = ylmo%dyn%now%cb_ref(i,j) * ngs%cf_south
            if (ylmo%bnd%basins(i,j) == 9.3_wp) ylmo%dyn%now%cb_ref(i,j) = ylmo%dyn%now%cb_ref(i,j) * ngs%cf_north
        end do
        end do

    end subroutine negis_update_cb_ref

    subroutine step_climate(dom, ts)
        ! Run climate on grid_clim and smb on grid_smb: geometry (z_srf/H_ice) from
        ! the hub, ocean/atmosphere forcing produced by snapclim, smb aggregated
        ! back to the Yelmo grid (conservative).
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        real(wp), allocatable :: z_srf_c(:,:), basins_c(:,:)
        character(len=256) :: gc, gn

        if (.not. dom%ctl%with_climate) return

        gc = trim(dom%ctl%grid_clim)
        gn = trim(dom%ctl%grid_name)

        ! snapclim snapshot on grid_clim, updated on the dt_clim cadence
        if (tstep_due(ts%time_elapsed, dom%ctl%dt_clim)) then
            call remap(dom, dom%topo%z_srf,   gn, z_srf_c,  gc, "bilin")
            call remap(dom, dom%topo%basins,  gn, basins_c, gc, "nn")
            call snapclim_update(dom%snp, z_srf=z_srf_c, time=ts%time, &
                                 domain=dom%ctl%domain, dx=dom%ctl%dx_clim, basins=basins_c)
        end if

        ! surface mass balance (smbpal or smb_simple), aggregated to the Yelmo grid
        call domain_update_smb(dom, ts)
    end subroutine step_climate

    subroutine domain_update_smb(dom, ts, init)
        ! Surface mass balance on grid_smb. Two methods: smbpal (default; monthly,
        ! needs tas/pr + geometry) or smb_simple (needs z_srf + sea-level
        ! temperature). Geometry comes from the hi-res hub, atmospheric forcing from
        ! snapclim (grid_clim). init=.true. runs the smbpal ITM equilibration before
        ! the first update. The result stays on grid_smb in the SMB model's fields
        ! (dom%smb%ann or dom%smbs); couple_smb_to_yelmo lands it on the Yelmo grid.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts
        logical, intent(in), optional    :: init

        real(wp), allocatable :: tas_s(:,:,:), pr_s(:,:,:), z_srf_s(:,:), H_ice_s(:,:)
        real(wp), allocatable :: tsl_s(:,:), Href_s(:,:)
        character(len=256) :: gc, gs, gn, gy
        logical :: is_init

        is_init = .false.
        if (present(init)) is_init = init

        gc = trim(dom%ctl%grid_clim)
        gs = trim(dom%ctl%grid_smb)
        gn = trim(dom%ctl%grid_name)
        gy = trim(dom%ctl%grid_yelmo)

        if (trim(dom%ctl%smb_method) == "smb_simple") then
            ! smb_simple: surface elevation + sea-level temperature, masked to the
            ! reference ice extent (refreshed each call in case H_ice_ref changed).
            call remap(dom, dom%topo%z_srf,          gn, z_srf_s, gs, "bilin")
            call remap(dom, dom%snp%now%tsl_ann,      gc, tsl_s,   gs, "bilin")
            call remap(dom, dom%yelmo%bnd%H_ice_ref,  gy, Href_s,  gs, "bilin")
            call smb_simple_set_mask(dom%smbs, Href_s)
            call smb_simple_update(dom%smbs, z_srf_s, tsl_s)
        else
            ! smbpal (monthly)
            call remap(dom, dom%snp%now%tas, gc, tas_s, gs, "bilin")
            call remap(dom, dom%snp%now%pr,  gc, pr_s,  gs, "bilin")
            call remap(dom, dom%topo%z_srf,  gn, z_srf_s, gs, "bilin")
            call remap(dom, dom%topo%H_ice,  gn, H_ice_s, gs, "bilin")
            if (is_init .and. trim(dom%smb%par%abl_method) == "itm") then
                call smbpal_update_monthly_equil(dom%smb, tas_s, pr_s, z_srf_s, H_ice_s, &
                        ts%time_rel, time_equil=100.0_wp)
            end if
            call smbpal_update_monthly(dom%smb, tas_s, pr_s, z_srf_s, H_ice_s, ts%time_rel)
        end if
    end subroutine domain_update_smb

    subroutine calc_glacial_smb(smb, lat2D, ta_ann, ta_ann_pd)
        ! Reduce (scale up toward zero) negative surface mass balance during
        ! glacial conditions, above a latitude limit. Ported verbatim from
        ! yelmox.f90; the glacial index is derived from the domain-mean cooling.
        real(wp), intent(inout) :: smb(:,:)
        real(wp), intent(in)    :: lat2D(:,:)
        real(wp), intent(in)    :: ta_ann(:,:)
        real(wp), intent(in)    :: ta_ann_pd(:,:)

        integer  :: i, j, nx, ny
        real(wp) :: t0, tnow, at
        real(wp), parameter :: dt_lgm  = -8.0_wp
        real(wp), parameter :: lat_lim = 55.0_wp
        real(wp), parameter :: fac_lim = 0.90_wp

        nx = size(smb,1)
        ny = size(smb,2)

        ! Quasi glacial-interglacial index (0: interglacial, 1: glacial)
        tnow = sum(ta_ann)    / real(nx*ny,wp)
        t0   = sum(ta_ann_pd) / real(nx*ny,wp)
        at = (tnow-t0)/dt_lgm
        if (at .lt. 0.0_wp) at = 0.0_wp
        if (at .gt. 1.0_wp) at = 1.0_wp

        do j = 1, ny
        do i = 1, nx
            if (smb(i,j) .lt. 0.0_wp .and. lat2D(i,j) .gt. lat_lim) then
                smb(i,j) = smb(i,j) - smb(i,j) * at * fac_lim
            end if
        end do
        end do
    end subroutine calc_glacial_smb

    subroutine refresh_htopo(dom)
        ! Refresh the hi-res geometry hub from the prognostic models (Yelmo grid
        ! -> hub grid, bilinear). The hub is then the geometry source for the
        ! coupling steps. Static masks (regions/basins) are not refreshed.
        type(ice_domain), intent(inout) :: dom

        call remap(dom, dom%yelmo%tpo%now%H_ice,  dom%ctl%grid_yelmo, &
                              dom%topo%H_ice,  dom%ctl%grid_name, "bilin")
        call remap(dom, dom%yelmo%bnd%z_bed,      dom%ctl%grid_yelmo, &
                              dom%topo%z_bed,  dom%ctl%grid_name, "bilin")
        call remap(dom, dom%yelmo%tpo%now%f_grnd, dom%ctl%grid_yelmo, &
                              dom%topo%f_grnd, dom%ctl%grid_name, "bilin")
        call remap(dom, dom%yelmo%bnd%z_sl,       dom%ctl%grid_yelmo, &
                              dom%topo%z_sl,   dom%ctl%grid_name, "bilin")
        call remap(dom, dom%yelmo%tpo%now%z_srf,  dom%ctl%grid_yelmo, &
                              dom%topo%z_srf,  dom%ctl%grid_name, "bilin")
    end subroutine refresh_htopo

    subroutine step_marine_shelf(dom, ts)
        ! Run marine_shelf on its own grid: geometry/masks from the hub, ocean
        ! forcing from snapclim. The outputs stay on grid_mshlf (in dom%mshlf%now);
        ! couple_marine_to_yelmo lands bmb_shlf / T_shlf on the Yelmo grid.
        type(ice_domain),  intent(inout) :: dom
        type(tstep_class), intent(in)    :: ts

        real(wp), allocatable :: H_ice_m(:,:), z_bed_m(:,:), f_grnd_m(:,:), z_sl_m(:,:)
        real(wp), allocatable :: regions_m(:,:), basins_m(:,:)
        real(wp), allocatable :: to_m(:,:,:), so_m(:,:,:), dto_m(:,:,:), dto_y(:,:,:)
        character(len=256) :: gm, gn, gc

        if (.not. dom%ctl%with_marine_shelf) return

        gm = trim(dom%ctl%grid_mshlf)
        gn = trim(dom%ctl%grid_name)
        gc = trim(dom%ctl%grid_clim)

        ! geometry + masks: hub -> mshlf grid
        call remap(dom, dom%topo%H_ice,   gn, H_ice_m,   gm, "bilin")
        call remap(dom, dom%topo%z_bed,   gn, z_bed_m,   gm, "bilin")
        call remap(dom, dom%topo%f_grnd,  gn, f_grnd_m,  gm, "bilin")
        call remap(dom, dom%topo%z_sl,    gn, z_sl_m,    gm, "bilin")
        call remap(dom, dom%topo%regions, gn, regions_m, gm, "nn")
        call remap(dom, dom%topo%basins,  gn, basins_m,  gm, "nn")

        ! ocean forcing (3D): snapclim (grid_clim) -> mshlf grid
        call remap(dom, dom%snp%now%to_ann, gc, to_m, gm, "bilin")
        call remap(dom, dom%snp%now%so_ann, gc, so_m, gm, "bilin")
        dto_y = dom%snp%now%to_ann - dom%snp%clim0%to_ann
        call remap(dom, dto_y, gc, dto_m, gm, "bilin")

        ! run marine_shelf on grid_mshlf
        call marshelf_update_shelf(dom%mshlf, H_ice_m, z_bed_m, f_grnd_m, basins_m, z_sl_m, &
                dom%ctl%dx_mshlf, dom%snp%now%depth, to_m, so_m, dto_ann=dto_m)
        call marshelf_update(dom%mshlf, H_ice_m, z_bed_m, f_grnd_m, regions_m, basins_m, &
                z_sl_m, dx=dom%ctl%dx_mshlf)
    end subroutine step_marine_shelf

    ! ----- output (climber-x convention: one file per module, on its own grid,
    !        named <module>_<grid>.nc; 1D timeseries as <module>_<grid>_ts.nc) ---

    function io_fname(outfldr, base, grid) result(fnm)
        character(len=*), intent(in) :: outfldr, base, grid
        character(len=512) :: fnm
        fnm = trim(outfldr)//trim(base)//"_"//trim(grid)//".nc"
    end function io_fname

    subroutine domain_write_init(dom, outfldr, time)
        ! Create the enabled per-module 2D output files (dims + static fields).
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr
        real(wp),         intent(in)    :: time

        if (dom%ctl%write_yelmo) &
            call yelmo_write_init(dom%yelmo, trim(io_fname(outfldr,"yelmo",dom%ctl%grid_yelmo)), &
                                  time_init=time, units="years")
        if (dom%ctl%write_htopo) &
            call htopo_write_init(dom%topo, trim(io_fname(outfldr,"htopo",dom%ctl%grid_name)), time_init=time)
        if (dom%ctl%write_isos) &
            call io_dims_init(trim(io_fname(outfldr,"isos",dom%ctl%grid_isos)),  dom%ctl%grid_isos,  time)
        if (dom%ctl%write_mshlf) &
            call io_dims_init(trim(io_fname(outfldr,"mshlf",dom%ctl%grid_mshlf)), dom%ctl%grid_mshlf, time)
        if (dom%ctl%write_smb) &
            call io_dims_init(trim(io_fname(outfldr,"smbpal",dom%ctl%grid_smb)),  dom%ctl%grid_smb,   time)
        if (dom%ctl%write_snap) &
            call io_dims_init(trim(io_fname(outfldr,"snap",dom%ctl%grid_clim)),   dom%ctl%grid_clim,  time)
    end subroutine domain_write_init

    subroutine domain_write_step(dom, outfldr, time)
        ! Append one time record to each enabled per-module 2D output file.
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr
        real(wp),         intent(in)    :: time

        ! Yelmo 2D fields: the yelmo_write_step default set + bmb_shlf, the
        ! coupled shelf melt (yelmo%bnd%bmb_shlf) aggregated onto the Yelmo grid.
        character(len=56), parameter :: yelmo_vars(23) = [ character(len=56) :: &
            "H_ice","z_srf","z_bed","mask_bed","uxy_b","uxy_s","uxy_bar", &
            "ux_bar","uy_bar","cb_ref","N_eff","beta","taub","taud","visc_bar", &
            "T_prime_b","hyd_W_til","mb_net","smb","bmb","cmb","z_sl","bmb_shlf" ]

        if (dom%ctl%write_yelmo) &
            call yelmo_write_step(dom%yelmo, trim(io_fname(outfldr,"yelmo",dom%ctl%grid_yelmo)), &
                                  time, nms=yelmo_vars, compare_pd=.FALSE.)
        if (dom%ctl%write_htopo) &
            call htopo_write_step(dom%topo, trim(io_fname(outfldr,"htopo",dom%ctl%grid_name)), time)
        if (dom%ctl%write_isos) &
            call isos_write_step(dom%isos, trim(io_fname(outfldr,"isos",dom%ctl%grid_isos)), time)
        if (dom%ctl%write_mshlf) &
            call mshlf_write_step(dom%mshlf, trim(io_fname(outfldr,"mshlf",dom%ctl%grid_mshlf)), time)
        if (dom%ctl%write_smb) &
            call smb_write_step(dom%smb, trim(io_fname(outfldr,"smbpal",dom%ctl%grid_smb)), time)
        if (dom%ctl%write_snap) &
            call snap_write_step(dom%snp, trim(io_fname(outfldr,"snap",dom%ctl%grid_clim)), time)
    end subroutine domain_write_step

    subroutine domain_write_1D(dom, outfldr, time, init)
        ! Write 1D timeseries: Yelmo regional aggregates + isostasy diagnostics.
        type(ice_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: outfldr
        real(wp),         intent(in)    :: time
        logical, intent(in), optional   :: init

        logical :: is_init
        character(len=512) :: fnm_isos

        is_init = .false.
        if (present(init)) is_init = init

        if (dom%ctl%write_yelmo) then
            if (is_init) then
                call yelmo_regions_write(dom%yelmo, time, init=.TRUE., units="years")
            else
                call yelmo_regions_write(dom%yelmo, time)
            end if
        end if

        if (dom%ctl%write_isos) then
            fnm_isos = trim(outfldr)//"isos_"//trim(dom%ctl%grid_isos)//"_ts.nc"
            if (is_init) call isos_write_1D_init(trim(fnm_isos), time)
            call isos_write_1D_step(dom%isos, trim(fnm_isos), time)
        end if
    end subroutine domain_write_1D

    ! --- private output helpers ---

    subroutine io_dims_init(filename, grid_name, time)
        ! Create a 2D output file with xc/yc (from the grid table) + time dims.
        character(len=*), intent(in) :: filename, grid_name
        real(wp),         intent(in) :: time
        type(grid_class) :: g
        call grid_cdo_read_desc(g, trim(grid_name), MAP_FLDR)
        call nc_create(filename)
        call nc_write_dim(filename, "xc", x=g%G%x, units="km")
        call nc_write_dim(filename, "yc", x=g%G%y, units="km")
        call nc_write_dim(filename, "time", x=time, dx=1.0_wp, nx=1, units="year", unlimited=.TRUE.)
    end subroutine io_dims_init

    subroutine io_var2D(filename, vnm, var, n, ncid, units, long_name)
        character(len=*), intent(in) :: filename, vnm, units, long_name
        real(wp),         intent(in) :: var(:,:)
        integer,          intent(in) :: n, ncid
        call nc_write(filename, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[size(var,1),size(var,2),1], ncid=ncid, &
                      units=units, long_name=long_name)
    end subroutine io_var2D

    subroutine io_ts(filename, vnm, val, n, ncid, units, long_name)
        character(len=*), intent(in) :: filename, vnm, units, long_name
        real(wp),         intent(in) :: val
        integer,          intent(in) :: n, ncid
        call nc_write(filename, vnm, val, dim1="time", start=[n], count=[1], ncid=ncid, &
                      units=units, long_name=long_name)
    end subroutine io_ts

    subroutine isos_write_step(isos, filename, time)
        type(isos_class), intent(in) :: isos
        character(len=*), intent(in) :: filename
        real(wp),         intent(in) :: time
        integer :: ncid, n
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call io_var2D(filename, "z_bed", isos%out%z_bed, n, ncid, "m", "Bedrock elevation")
        call io_var2D(filename, "z_ss",  isos%out%z_ss,  n, ncid, "m", "Sea-surface height")
        call io_var2D(filename, "w",     isos%out%w,     n, ncid, "m", "Viscous displacement")
        call io_var2D(filename, "we",    isos%out%we,    n, ncid, "m", "Elastic displacement")
        call nc_close(ncid)
    end subroutine isos_write_step

    subroutine mshlf_write_step(mshlf, filename, time)
        type(marshelf_class), intent(in) :: mshlf
        character(len=*),     intent(in) :: filename
        real(wp),             intent(in) :: time
        integer :: ncid, n
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call io_var2D(filename, "bmb_shlf", mshlf%now%bmb_shlf, n, ncid, "m/yr", "Shelf basal mass balance")
        call io_var2D(filename, "T_shlf",   mshlf%now%T_shlf,   n, ncid, "K", "Shelf temperature")
        call io_var2D(filename, "tf_shlf",  mshlf%now%tf_shlf,  n, ncid, "K", "Thermal forcing")
        call nc_close(ncid)
    end subroutine mshlf_write_step

    subroutine smb_write_step(smb, filename, time)
        type(smbpal_class), intent(in) :: smb
        character(len=*),   intent(in) :: filename
        real(wp),           intent(in) :: time
        integer :: ncid, n
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call io_var2D(filename, "smb",  smb%ann%smb,  n, ncid, "m ie/yr", "Surface mass balance")
        call io_var2D(filename, "tsrf", smb%ann%tsrf, n, ncid, "K", "Surface temperature")
        call nc_close(ncid)
    end subroutine smb_write_step

    subroutine snap_write_step(snp, filename, time)
        type(snapclim_class), intent(in) :: snp
        character(len=*),     intent(in) :: filename
        real(wp),             intent(in) :: time
        integer :: ncid, n
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call io_var2D(filename, "t2m_ann", snp%now%ta_ann, n, ncid, "K", "Annual mean air temperature")
        call io_var2D(filename, "pr_ann",  snp%now%pr_ann, n, ncid, "mm/a", "Annual mean precipitation")
        call nc_close(ncid)
    end subroutine snap_write_step

    subroutine isos_write_1D_init(filename, time)
        character(len=*), intent(in) :: filename
        real(wp),         intent(in) :: time
        call nc_create(filename)
        call nc_write_dim(filename, "time", x=time, dx=1.0_wp, nx=1, units="year", unlimited=.TRUE.)
    end subroutine isos_write_1D_init

    subroutine isos_write_1D_step(isos, filename, time)
        type(isos_class), intent(in) :: isos
        character(len=*), intent(in) :: filename
        real(wp),         intent(in) :: time
        integer :: ncid, n, np
        call nc_open(filename, ncid, writable=.TRUE.)
        n = nc_time_index(filename, "time", time, ncid)
        np = size(isos%out%z_bed)
        call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call io_ts(filename, "bsl",        isos%now%bsl,                       n, ncid, "m", "Barystatic sea level")
        call io_ts(filename, "z_bed_mean", sum(isos%out%z_bed)/real(np,wp),    n, ncid, "m", "Mean bedrock elevation")
        call io_ts(filename, "z_bed_min",  minval(isos%out%z_bed),             n, ncid, "m", "Min bedrock elevation")
        call io_ts(filename, "z_bed_max",  maxval(isos%out%z_bed),             n, ncid, "m", "Max bedrock elevation")
        call io_ts(filename, "w_mean",     sum(isos%out%w)/real(np,wp),        n, ncid, "m", "Mean viscous displacement")
        call io_ts(filename, "we_mean",    sum(isos%out%we)/real(np,wp),       n, ncid, "m", "Mean elastic displacement")
        call nc_close(ncid)
    end subroutine isos_write_1D_step

    ! ----- remap: identity-copy when src == dst, else via the coupler ---

    subroutine remap_2D(dom, var_src, src, var_dst, dst, method)
        type(ice_domain),      intent(inout) :: dom
        real(wp),              intent(in)    :: var_src(:,:)
        character(len=*),      intent(in)    :: src, dst, method
        real(wp), allocatable, intent(inout) :: var_dst(:,:)

        if (trim(src) == trim(dst)) then
            if (allocated(var_dst)) then
                if (size(var_dst,1) /= size(var_src,1) .or. &
                    size(var_dst,2) /= size(var_src,2)) deallocate(var_dst)
            end if
            if (.not. allocated(var_dst)) allocate(var_dst(size(var_src,1), size(var_src,2)))
            var_dst = var_src
        else
            call cpl_remap(dom%cpl, var_src, src, var_dst, dst, method=method)
        end if
    end subroutine remap_2D

    subroutine remap_3D(dom, var_src, src, var_dst, dst, method)
        type(ice_domain),      intent(inout) :: dom
        real(wp),              intent(in)    :: var_src(:,:,:)
        character(len=*),      intent(in)    :: src, dst, method
        real(wp), allocatable, intent(inout) :: var_dst(:,:,:)

        if (trim(src) == trim(dst)) then
            if (allocated(var_dst)) then
                if (size(var_dst,1) /= size(var_src,1) .or. &
                    size(var_dst,2) /= size(var_src,2) .or. &
                    size(var_dst,3) /= size(var_src,3)) deallocate(var_dst)
            end if
            if (.not. allocated(var_dst)) &
                allocate(var_dst(size(var_src,1), size(var_src,2), size(var_src,3)))
            var_dst = var_src
        else
            call cpl_remap(dom%cpl, var_src, src, var_dst, dst, method=method)
        end if
    end subroutine remap_3D

end module yelmox_domain
