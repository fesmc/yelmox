module kryos
    ! One region's full model state, the kryos_domain: Yelmo plus its component
    ! models (isostasy, climate, smb, marine shelf, sediments, geothermal), the
    ! hi-res topography hub (dom%topo) and the coupler holding the grids and
    ! remap maps (dom%cpl).
    !
    ! The hub is the geometry source of truth: it is refreshed from the
    ! prognostic models each step (refresh_hub), and the coupling steps remap
    ! fields to/from it through the coupler. Each component runs on its own grid,
    ! named in domain_ctl; a component's grid is just a string, so e.g.
    ! marine_shelf can run on the hub grid or on the Yelmo grid simply by
    ! changing grid_mshlf.
    !
    ! This module defines the domain and its configuration, initializes it
    ! (domain_init) and provides remap, the domain-level grid transfer. The rest
    ! of the framework builds on it:
    !   kryos_regions   region-specific masks, sub-regions and physics
    !   kryos_coupling  the per-step coupling primitives (step_*, couple_*)
    !   kryos_startup   initial state (cold start) and restart bundles
    !   kryos_output    per-module 2D/1D output
    !   kryos_forcing   driver-owned transient forcing (tsgen)

    use nml,          only : nml_read, nml_replace
    use coords,       only : grid_class, grid_cdo_read_desc
    use yelmo,        only : yelmo_class, wp, yelmo_init, yelmo_init_grid, ytopo_input_class
    use yelmo_defs,   only : MASK_ICE_NONE, MASK_ICE_DYNAMIC
    use marine_shelf, only : marshelf_class, marshelf_init
    use fastisostasy, only : isos_class, isos_init
    use climate_out,    only : climate_out_class
    use yelmox_climate, only : yelmox_climate_class, climate_init
    use smbpal,       only : smbpal_class, smbpal_init
    use smb_simple_m, only : smb_simple_class, smb_simple_init, smb_simple_set_mask
    use ice_optimization, only : ice_opt_params, optimize_par_load
    use sediments,    only : sediments_class, sediments_init
    use geothermal,   only : geothermal_class, geothermal_init
    use htopo,        only : htopo_class, htopo_init, htopo_ice_allowed, htopo_relax_tau
    use coupler,      only : coupler_class, coupler_init, coupler_prime, cpl_remap => remap

    implicit none
    private

    ! Folder holding the grid descriptions (grid_<name>.txt) and cached maps.
    character(len=*), parameter :: MAP_FLDR = "maps"

    type domain_ctl
        ! Parameter file (kept for sub-steps that reload from it, e.g. LGM startup).
        character(len=512) :: path_par = ""
        ! Shared timeline (read from the driver's timeline group, e.g. [ctrl]).
        character(len=56) :: tstep_method = "const"
        real(wp) :: dtt         = 10.0_wp
        ! Cadences + methods ([coupling]).
        real(wp) :: dt_clim     = 10.0_wp   ! [yr] snapclim snapshot update frequency
        character(len=56) :: equil_method = "none"
        character(len=56) :: smb_method   = "smbpal"

        ! Cold-start ice state ([coupling]): init_marine_H first, then init_method.
        character(len=56)  :: init_method     = "equil"   ! none | equil | recon | recon_ref
        logical            :: init_marine_H   = .false.   ! impose LGM-like marine ice
        real(wp)           :: init_equil_time = 10.0_wp   ! [yr] equil: equilibration time
        character(len=512) :: recon_path = ""             ! recon*: ice reconstruction file
        character(len=56)  :: recon_var  = ""             ! recon*: its ice-thickness variable
        real(wp)           :: recon_codes(20)             ! recon: regions where its ice is imposed
        integer            :: n_recon_codes = 0

        ! Optional physics switches ([coupling]).
        logical :: scale_glacial_smb       = .false.   ! reduce negative glacial smb ([glacial_smb] group)
        logical :: lim_pd_ice              = .false.   ! extra melt outside PD ice extent
        logical :: use_negis               = .false.   ! NEGIS cb_ref modification ([negis] group)

        ! Which components are active in this domain's coupling sequence.
        logical :: with_ice_sheet    = .true.
        logical :: with_isostasy     = .true.
        logical :: with_marine_shelf = .true.
        logical :: with_climate      = .true.

        ! Domain name + grid of every component ([domain]; the source of truth
        ! for remap keys). A blank component grid takes its default.
        character(len=256) :: domain     = ""   ! e.g. "Antarctica"
        character(len=256) :: grid_hub   = ""   ! hi-res geometry hub, highest res
        character(len=256) :: grid_ice   = ""   ! Yelmo grid (default = grid_hub)
        character(len=256) :: grid_mshlf = ""   ! marine-shelf grid (default = grid_hub)
        real(wp) :: dx_mshlf = 0.0_wp           ! marine-shelf grid spacing (Yelmo dx units)
        character(len=256) :: grid_isos = ""    ! isostasy grid (default = grid_ice)
        real(wp) :: dx_isos = 0.0_wp            ! isostasy grid spacing in x (Yelmo dx units)
        real(wp) :: dy_isos = 0.0_wp            ! isostasy grid spacing in y (Yelmo dy units)
        ! grid_clim sets the grid of BOTH the reference climatology (often from a
        ! high-resolution regional model) and the transient forcing (often from
        ! a coarser climate model). A coarse grid_clim matches the forcing but
        ! loses detail of the high-res reference. Until the two get separate
        ! grids, set grid_clim to the highest-resolution climate input.
        character(len=256) :: grid_clim = ""    ! climate grid (default = grid_ice)
        real(wp) :: dx_clim = 0.0_wp            ! climate grid spacing (Yelmo dx units)
        character(len=256) :: grid_smb = ""     ! smb grid (default = grid_clim)

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
        real(wp) :: basin_centre = 9.1_wp   ! basin codes of the NEGIS parts
        real(wp) :: basin_south  = 9.2_wp
        real(wp) :: basin_north  = 9.3_wp
    end type negis_params

    ! Glacial smb scaling parameters: negative smb above lat_lim is reduced by
    ! up to fac_lim, with a glacial index from the domain-mean cooling (dt_lgm = full glacial).
    type glacial_smb_params
        real(wp) :: dt_lgm  = -8.0_wp    ! [K] domain-mean cooling of a full glacial
        real(wp) :: lat_lim = 55.0_wp    ! [deg] latitude above which smb is scaled
        real(wp) :: fac_lim = 0.9_wp     ! [1] maximum reduction of negative smb
    end type glacial_smb_params

    type kryos_domain
        type(yelmo_class)      :: yelmo
        type(marshelf_class)   :: mshlf
        type(isos_class)       :: isos
        type(yelmox_climate_class) :: cl    ! climate backend (snapclim | snapesm)
        type(climate_out_class)    :: clim  ! backend-agnostic climate output (now/ref)
        type(smbpal_class)     :: smb
        type(smb_simple_class) :: smbs    ! alternative SMB (smb_method="smb_simple")
        type(sediments_class)  :: sed
        type(geothermal_class) :: gthrm
        type(htopo_class)      :: topo    ! hi-res geometry reference hub
        type(coupler_class)    :: cpl     ! this region's grid resolution + map cache
        type(ice_opt_params)   :: opt     ! basal-friction / thermal-forcing optimization
        type(negis_params)     :: ngs     ! NEGIS cb_ref modification
        type(glacial_smb_params) :: gsmb  ! glacial smb scaling
        type(domain_ctl)       :: ctl
    end type kryos_domain

    public :: MAP_FLDR
    public :: domain_ctl, negis_params, glacial_smb_params, kryos_domain
    public :: domain_init
    public :: cadence_due
    ! remap is used by every coupling step, and by flavor drivers (e.g. the ESM
    ! driver) to move fields between the hub/Yelmo grids and their own grid.
    public :: remap, remap_method_smooth

    ! Domain-level remap: identity-copy when src == dst, else remap via the coupler.
    interface remap
        module procedure remap_2D, remap_3D
    end interface remap

contains

    function cadence_due(time, dt) result(due)
        ! Cadence predicate: true when `time` falls on the dt grid (0.01-yr
        ! precision). dt <= 0 disables the cadence (never due).
        real(wp), intent(in) :: time, dt
        logical :: due
        due = .false.
        if (dt > 0.0_wp) due = (mod(nint(time*100), nint(dt*100)) == 0)
    end function cadence_due

    subroutine domain_init(dom, path_par, time, group_suffix, init_climate, timeline_group)
        ! Initialize all sub-models of one domain, load the hi-res reference hub,
        ! prime the Yelmo<->hub maps, and place marine_shelf on its configured grid.
        ! The barystatic sea level (bsl) is NOT a domain sub-model: it is a shared,
        ! driver-owned object (one per run, common to every domain), so it is
        ! initialized by the driver and passed into the isostasy steps.
        !
        ! The domain definition ([domain]: name, component grids, hub topography
        ! and masks) is read first and the hub loaded. Yelmo is then populated
        ! from the domain like the other components: its grid, and its initial
        ! and present-day topography remapped from the hub.
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
        ! tstep_init -- from which the domain reads tstep_method/dtt itself.
        type(kryos_domain), intent(inout) :: dom
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
        real(wp), allocatable :: regions_y(:,:), basins_y(:,:)
        integer,  allocatable :: mask_ice_y(:,:)
        real(wp), allocatable :: xs(:), ys(:), lats_s(:,:), Href_s(:,:)
        type(ytopo_input_class) :: topo_y

        sfx = ""
        if (present(group_suffix)) sfx = trim(group_suffix)

        do_climate = .TRUE.
        if (present(init_climate)) do_climate = init_climate

        tgroup = "ctrl"
        if (present(timeline_group)) tgroup = trim(timeline_group)

        ! --- domain definition + run control ---
        call domain_ctl_load(dom%ctl, path_par, trim(sfx), trim(tgroup))
        domain = trim(dom%ctl%domain)

        ! --- hi-res geometry hub (topography + masks from [domain]) + coupler ---
        call htopo_init(dom%topo, path_par, "domain"//trim(sfx), domain, dom%ctl%grid_hub)

        ! Grids resolve from maps/grid_<name>.txt; prime the Yelmo<->hub maps.
        call coupler_init(dom%cpl)
        call coupler_prime(dom%cpl, dom%ctl%grid_ice, dom%ctl%grid_hub, "bilin")  ! Yelmo -> hub
        call coupler_prime(dom%cpl, dom%ctl%grid_hub, dom%ctl%grid_ice, "con")    ! hub -> Yelmo
        call coupler_prime(dom%cpl, dom%ctl%grid_hub, dom%ctl%grid_ice, "nn")     ! hub -> Yelmo (masks)

        ! --- ice sheet on grid_ice, with the hub's topography and masks ---
        ! The hub topography is both Yelmo's initial state and its present-day
        ! reference (H_ice_ref, z_bed_ref, optimization target). The code masks
        ! come from the hub too, and the domain says where ice is allowed
        ! (Yelmo's mask_border then sets the border).
        call grid_cdo_read_desc(grid_y, trim(dom%ctl%grid_ice), MAP_FLDR)
        call yelmo_init_grid(dom%yelmo%grd, grid_y)

        call remap(dom, dom%topo%z_bed,    dom%ctl%grid_hub, topo_y%z_bed,    dom%ctl%grid_ice, "con")
        call remap(dom, dom%topo%H_ice,    dom%ctl%grid_hub, topo_y%H_ice,    dom%ctl%grid_ice, "con")
        call remap(dom, dom%topo%z_srf,    dom%ctl%grid_hub, topo_y%z_srf,    dom%ctl%grid_ice, "con")
        call remap(dom, dom%topo%z_bed_sd, dom%ctl%grid_hub, topo_y%z_bed_sd, dom%ctl%grid_ice, "con")

        call remap(dom, dom%topo%regions,  dom%ctl%grid_hub, regions_y,       dom%ctl%grid_ice, "nn")
        call remap(dom, dom%topo%basins,   dom%ctl%grid_hub, basins_y,        dom%ctl%grid_ice, "nn")
        allocate(mask_ice_y(size(regions_y,1), size(regions_y,2)))
        mask_ice_y = MASK_ICE_NONE
        where (htopo_ice_allowed(dom%topo%par, regions_y)) mask_ice_y = MASK_ICE_DYNAMIC

        call yelmo_init(dom%yelmo, filename=path_par, grid_def="none", time=time, &
                        domain=domain, grid_name=dom%ctl%grid_ice, &
                        group="yelmo"//trim(sfx), regions=regions_y, basins=basins_y, &
                        mask_ice=mask_ice_y, topo_init=topo_y, topo_pd=topo_y)

        ! Where the ice relaxes to the reference (ytopo.topo_rel = -1).
        dom%yelmo%bnd%tau_relax = htopo_relax_tau(dom%topo%par, regions_y)

        ! --- external forcing models (climate/smb/isostasy on the Yelmo grid) ---
        ! Isostasy on its configured grid (grid_isos).
        call grid_cdo_read_desc(grid_i, trim(dom%ctl%grid_isos),  MAP_FLDR)
        nx_i = grid_i%G%nx
        ny_i = grid_i%G%ny
        ! Grid spacing in Yelmo units, scaled by the resolution ratio (per axis).
        dom%ctl%dx_isos = dom%yelmo%grd%G%dx * (grid_i%G%dx / grid_y%G%dx)
        dom%ctl%dy_isos = dom%yelmo%grd%G%dy * (grid_i%G%dy / grid_y%G%dy)
        call isos_init(dom%isos, path_par, "isos"//trim(sfx), nx_i, ny_i, &
                       dom%ctl%dx_isos, dom%ctl%dy_isos, cnst=dom%yelmo%bnd%cnst)

        call sediments_init(dom%sed, path_par, dom%yelmo%grd%G%nx, dom%yelmo%grd%G%ny, &
                            domain, dom%ctl%grid_ice, group="sed"//trim(sfx))
        dom%yelmo%bnd%H_sed = dom%sed%now%H

        call geothermal_init(dom%gthrm, path_par, dom%yelmo%grd%G%nx, dom%yelmo%grd%G%ny, &
                             domain, dom%ctl%grid_ice, group="ghf"//trim(sfx))
        dom%yelmo%bnd%Q_geo = dom%gthrm%now%ghf

        ! --- climate on its configured grid (grid_clim) ---
        ! snapclim reads grid-specific input data, so grid_clim must be a grid whose
        ! forcing files exist (the Yelmo grid for the standard setup).
        call grid_cdo_read_desc(grid_c, trim(dom%ctl%grid_clim), MAP_FLDR)
        nx_c = grid_c%G%nx
        ny_c = grid_c%G%ny
        dom%ctl%dx_clim = dom%yelmo%grd%G%dx * (grid_c%G%dx / grid_y%G%dx)
        if (do_climate) then
            call remap(dom, dom%topo%basins, dom%ctl%grid_hub, basins_c, dom%ctl%grid_clim, "nn")
            call climate_init(dom%cl, path_par, domain, trim(dom%ctl%grid_clim), &
                              nx_c, ny_c, time, basins_c, group="snap"//trim(sfx))
        end if

        ! --- smb on its configured grid (grid_smb) ---
        ! smbpal reads no grid-specific data; only lats (insolation) is physical.
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
                                 group="smb_simple"//trim(sfx), units="m", &
                                 cnst=dom%yelmo%bnd%cnst)
            call remap(dom, dom%yelmo%bnd%H_ice_ref, dom%ctl%grid_ice, &
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
        call remap(dom, dom%topo%regions, dom%ctl%grid_hub, regions_m, dom%ctl%grid_mshlf, "nn")
        call remap(dom, dom%topo%basins,  dom%ctl%grid_hub, basins_m,  dom%ctl%grid_mshlf, "nn")

        call marshelf_init(dom%mshlf, path_par, "marine_shelf"//trim(sfx), nx_m, ny_m, &
                           domain, trim(dom%ctl%grid_mshlf), regions_m, basins_m, &
                           cnst=dom%yelmo%bnd%cnst)

        ! Optimization state (basal friction + thermal forcing); no-op unless
        ! equil_method == "opt". Must follow yelmo_init (grid + till params known).
        call domain_opt_init(dom, path_par, trim(sfx))

        ! NEGIS cb_ref modification and glacial smb scaling: load their groups
        ! when enabled, so they cannot silently run with default parameters.
        if (dom%ctl%use_negis)         call negis_par_load(dom%ngs, path_par, trim(sfx))
        if (dom%ctl%scale_glacial_smb) call glacial_smb_par_load(dom%gsmb, path_par, trim(sfx))

    end subroutine domain_init

    subroutine negis_par_load(ngs, path_par, suffix)
        ! Load the NEGIS cb_ref scaling parameters ([negis<suffix>]). Only read
        ! when [coupling] use_negis is set.
        type(negis_params), intent(inout) :: ngs
        character(len=*),   intent(in)    :: path_par
        character(len=*),   intent(in)    :: suffix

        ngs%use_negis_par = .true.
        call nml_read(path_par, "negis"//trim(suffix), "cf_0",      ngs%cf_0)
        call nml_read(path_par, "negis"//trim(suffix), "cf_1",      ngs%cf_1)
        call nml_read(path_par, "negis"//trim(suffix), "cf_centre", ngs%cf_centre)
        call nml_read(path_par, "negis"//trim(suffix), "cf_north",  ngs%cf_north)
        call nml_read(path_par, "negis"//trim(suffix), "cf_south",  ngs%cf_south)
        call nml_read(path_par, "negis"//trim(suffix), "basin_centre", ngs%basin_centre)
        call nml_read(path_par, "negis"//trim(suffix), "basin_south",  ngs%basin_south)
        call nml_read(path_par, "negis"//trim(suffix), "basin_north",  ngs%basin_north)
    end subroutine negis_par_load

    subroutine glacial_smb_par_load(gsmb, path_par, suffix)
        ! Load the glacial smb scaling parameters ([glacial_smb<suffix>]). Only
        ! read when [coupling] scale_glacial_smb is set.
        type(glacial_smb_params), intent(inout) :: gsmb
        character(len=*),         intent(in)    :: path_par
        character(len=*),         intent(in)    :: suffix

        call nml_read(path_par, "glacial_smb"//trim(suffix), "dt_lgm",  gsmb%dt_lgm)
        call nml_read(path_par, "glacial_smb"//trim(suffix), "lat_lim", gsmb%lat_lim)
        call nml_read(path_par, "glacial_smb"//trim(suffix), "fac_lim", gsmb%fac_lim)
    end subroutine glacial_smb_par_load

    subroutine domain_opt_init(dom, path_par, suffix)
        ! Load optimization parameters and prepare Yelmo for external cb_ref:
        ! allocate/seed the friction bounds (cf_min/cf_max) on the Yelmo grid and
        ! switch till_method to external (-1) so yelmo_update uses the optimized
        ! cb_ref. The initial cb_ref guess (cold start) is set in domain_init_state
        ! after yelmo_init_state; on restart cb_ref is restored from the bundle.
        ! No-op unless equil_method == "opt".
        type(kryos_domain), intent(inout) :: dom
        character(len=*), intent(in)    :: path_par
        character(len=*), intent(in)    :: suffix

        integer :: nx, ny

        if (trim(dom%ctl%equil_method) /= "opt") return

        dom%opt%tf_basins = 0
        call optimize_par_load(dom%opt, path_par, "opt"//trim(suffix))

        if (dom%opt%cf_init <= 0.0_wp) then
            write(*,*) "domain_opt_init:: error: opt"//trim(suffix)//".cf_init must be > 0 &
                       &(the cold-start cb_ref guess); got ", dom%opt%cf_init
            stop 1
        end if

        nx = dom%yelmo%grd%G%nx
        ny = dom%yelmo%grd%G%ny
        allocate(dom%opt%cf_min(nx, ny), dom%opt%cf_max(nx, ny))
        dom%opt%cf_min = dom%yelmo%dyn%par%till_cf_min
        dom%opt%cf_max = dom%yelmo%dyn%par%till_cf_ref

        ! cb_ref is set externally by the optimization from here on.
        dom%yelmo%dyn%par%till_method = -1
    end subroutine domain_opt_init

    subroutine domain_ctl_load(ctl, path_par, suffix, timeline_group)
        ! Load this domain's definition + coupling + output config. All groups
        ! carry an optional domain suffix (e.g. "_north"), so several domains can
        ! coexist in one parameter file without group.name collisions (matters
        ! for runme -p). The shared timeline is driver-owned (tstep_init); the
        ! values the domain logic needs (tstep_method, dtt) are read from the
        ! same timeline_group here, so nothing is injected after init.
        type(domain_ctl), intent(inout) :: ctl
        character(len=*), intent(in)    :: path_par
        character(len=*), intent(in)    :: suffix
        character(len=*), intent(in)    :: timeline_group

        character(len=256) :: gd, gc, go

        ctl%path_par = trim(path_par)

        ! Shared timeline values used by the domain logic (dt_clim cadence,
        ! optimization dt, domain-specific startup).
        call nml_read(path_par, timeline_group, "tstep_method", ctl%tstep_method)
        call nml_read(path_par, timeline_group, "dtt",          ctl%dtt)

        ! Domain definition ([domain<suffix>]): name and the grid of every
        ! component. The hub's topography and masks are read by htopo_init.
        gd = "domain"//trim(suffix)
        call nml_read(path_par, gd, "name",     ctl%domain)
        call nml_read(path_par, gd, "grid_hub", ctl%grid_hub)
        ctl%grid_ice   = ""
        ctl%grid_isos  = ""
        ctl%grid_clim  = ""
        ctl%grid_smb   = ""
        ctl%grid_mshlf = ""
        call nml_read(path_par, gd, "grid_ice",   ctl%grid_ice)
        call nml_read(path_par, gd, "grid_isos",  ctl%grid_isos)
        call nml_read(path_par, gd, "grid_clim",  ctl%grid_clim)
        call nml_read(path_par, gd, "grid_smb",   ctl%grid_smb)
        call nml_read(path_par, gd, "grid_mshlf", ctl%grid_mshlf)

        ! A blank component grid takes its default.
        if (len_trim(ctl%grid_ice)   == 0) ctl%grid_ice   = trim(ctl%grid_hub)
        if (len_trim(ctl%grid_isos)  == 0) ctl%grid_isos  = trim(ctl%grid_ice)
        if (len_trim(ctl%grid_clim)  == 0) ctl%grid_clim  = trim(ctl%grid_ice)
        if (len_trim(ctl%grid_smb)   == 0) ctl%grid_smb   = trim(ctl%grid_clim)
        if (len_trim(ctl%grid_mshlf) == 0) ctl%grid_mshlf = trim(ctl%grid_hub)

        ! Coupling ([coupling<suffix>]): active components, methods, restart
        ! bundle. The restart cadence is the driver's [tm_rst] timeout.
        gc = "coupling"//trim(suffix)
        call nml_read(path_par, gc, "with_ice_sheet",    ctl%with_ice_sheet)
        call nml_read(path_par, gc, "with_isostasy",     ctl%with_isostasy)
        call nml_read(path_par, gc, "with_climate",      ctl%with_climate)
        call nml_read(path_par, gc, "with_marine_shelf", ctl%with_marine_shelf)
        call nml_read(path_par, gc, "equil_method",   ctl%equil_method)
        ctl%smb_method = "smbpal"
        call nml_read(path_par, gc, "smb_method",     ctl%smb_method)
        call nml_read(path_par, gc, "dt_clim",        ctl%dt_clim)

        ! Optional physics switches. use_negis and scale_glacial_smb
        ! additionally require a [negis<suffix>] / [glacial_smb<suffix>] group.
        call nml_read(path_par, gc, "scale_glacial_smb",       ctl%scale_glacial_smb)
        call nml_read(path_par, gc, "lim_pd_ice",              ctl%lim_pd_ice)
        call nml_read(path_par, gc, "use_negis",               ctl%use_negis)

        ! Cold-start ice state. The keys of the selected init_method are required.
        call nml_read(path_par, gc, "init_method",   ctl%init_method)
        call nml_read(path_par, gc, "init_marine_H", ctl%init_marine_H)
        select case(trim(ctl%init_method))
            case("none")
            case("equil")
                call nml_read(path_par, gc, "init_equil_time", ctl%init_equil_time)
            case("recon", "recon_ref")
                call nml_read(path_par, gc, "recon_path", ctl%recon_path)
                call nml_read(path_par, gc, "recon_var",  ctl%recon_var)
                ! {domain}/{grid_name} resolve to the domain name and grid_ice.
                call nml_replace(ctl%recon_path, "{domain}",    trim(ctl%domain))
                call nml_replace(ctl%recon_path, "{grid_name}", trim(ctl%grid_ice))
                if (trim(ctl%init_method) == "recon") then
                    ctl%recon_codes = -9999.0_wp
                    call nml_read(path_par, gc, "recon_codes", ctl%recon_codes)
                    ctl%n_recon_codes = count(ctl%recon_codes /= -9999.0_wp)
                    if (ctl%n_recon_codes == 0) then
                        write(*,*) "domain_ctl_load:: error: "//trim(gc)//".init_method = recon needs recon_codes."
                        stop 1
                    end if
                end if
            case default
                write(*,*) "domain_ctl_load:: error: "//trim(gc)//".init_method must be none, equil, &
                           &recon or recon_ref; got "//trim(ctl%init_method)
                stop 1
        end select

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

    ! ----- remap: identity-copy when src == dst, else via the coupler ---

    function remap_method_smooth(dx_src, dx_dst) result(method)
        ! Remap method for a field crossing the Yelmo/isostasy boundary: bilinear
        ! to refine (coarse -> fine), conservative to coarsen (fine -> coarse).
        ! Coarsening an ice load or a displacement by averaging over the target
        ! cell is what keeps the driving mass, rather than sampling one point of
        ! it; refining a long-wavelength field wants the smooth interpolant, not
        ! piecewise-constant blocks. Equal spacings never reach a real remap --
        ! remap_2D short-circuits to a copy.
        real(wp), intent(in) :: dx_src, dx_dst
        character(len=32) :: method

        if (dx_dst < dx_src) then
            method = "bilin"
        else
            method = "con"
        end if
    end function remap_method_smooth

    subroutine remap_2D(dom, var_src, src, var_dst, dst, method)
        type(kryos_domain),      intent(inout) :: dom
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
        type(kryos_domain),      intent(inout) :: dom
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

end module kryos
