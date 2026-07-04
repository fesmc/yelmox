# Multigrid coupling (`yelmox`)

Design doc for a multigrid rewrite of the yelmox driver. Status: **implemented**
(`yelmox.f90` + `libs/yelmox_domain.f90`); single-domain parity with
`yelmox.f90` validated, multi-domain (bipolar) runs, optimization + smb_simple +
domain-specific startups ported. The bipolar driver (`yelmox_bipolar.f90`) also
carries the full `yelmox_bipolar` ocean coupling: a shared barystatic sea level
and a shared Ocean Box Model exchanging freshwater flux / ocean temperature
between hemispheres (`yelmox_osm/obm_coupling.f90`). A third driver (`yelmox_esm.f90`)
swaps snapclim for ESM climatic forcing (`libs/esm.f90`), running esm as a
first-class component on its own grid (`grid_clim`) and remapping its outputs to
the consumer grids, just like snapclim (see Drivers). Remaining: FastIsostasy
hi-res output (below).

## Motivation

Today every helper module (`marine_shelf`, `snapclim`, `smbpal`, FastIsostasy)
is initialized on, and operates on, the Yelmo grid. Fields are passed between
modules as bare `(nx,ny)` 2D arrays shaped to that single grid. There is no
inter-module remapping in the main loop — remapping (`map_class` from the
`coords` library in fesm-utils) is used only for restart I/O.

We want each module to be able to work on its own grid/resolution and be
remapped to Yelmo (or any other module) at the moment of coupling. Example
workflow for Antarctica with Yelmo on `ANT-16KM` and a hi-res topo/working grid
`ANT-2KM`:

- `H_ice` (Yelmo `ANT-16KM`) → bilin → `ANT-2KM` for the marine-shelf calc.
- `z_bed`/`z_ss` from FastIsostasy → conservatively aggregated to Yelmo's grid,
  set as BCs.
- marine_shelf computed at `ANT-2KM`; `bmb_shlf`/`T_shlf` → conservatively
  aggregated → `ANT-16KM` → passed to Yelmo as forcing.

The bookkeeping of many remapping steps is the thing to tame.

## Key observations from the current code

- **Helper modules are already grid-agnostic.** `marshelf_update`,
  `isos_update`, `snapclim_update`, `smbpal_update` all take/return plain
  `(nx,ny)` arrays. They are tied to the Yelmo grid only because the *caller*
  feeds them Yelmo-shaped arrays. Moving a module onto another grid is a
  caller-side change, not a module change.
- **`map_class` self-identifies and self-caches.** It stores `name1`, `name2`,
  `method`, and `map_init` loads from / saves to a `maps/` folder on disk. A map
  object already carries its own identity key; disk caching across runs is free.
- **`method` is baked into the map** at `map_init` time (`"con"` = area weights,
  `"bilin"` = neighbor weights). So `method` is part of the map cache key.
  **`stat`** (`mean`/`count`/`stdev`) is applied at `map_field` time and does not
  change the map — it is a per-call pass-through, not a key.
- **`map_field` is generic over kind.** The interface already covers `dp`, `sp`
  (`map_field_grid_grid_sp`) and `int` (`map_field_grid_grid_int`); the sp/int
  variants accumulate in `dp` internally. So `remap` never converts precision
  outside — it calls the generic `map_field` and lets it dispatch by kind.

## Conceptual model

Grids are **nodes**; maps are **directed, method-typed edges**. `16KM→2KM bilin`
(downscale) and `2KM→16KM con` (aggregate) are two distinct objects even between
the same pair. Worst case is `N·(N−1)·methods` edges, but a real run uses a
handful. Maps are built **lazily** (only edges actually used) and each is built
**once**, stored **once**, and **shared** — the only hard requirement, because a
hi-res target map is both expensive to build and large in memory.

## Layering

```
fesm-utils/utils/src/coords/coupler.f90    coupler_class + remap
        │  use coords  (grid_class, grid_init, map_class, map_init, map_field)
        ▼
yelmox/libs/yelmox_domain.f90              ice_domain + step_* + yelmox_step
        │  use coupler, yelmo, marine_shelf, fastisostasy, snapclim, smbpal, bsl
        ▼
yelmox/yelmox.f90                        thin driver: init → time loop → I/O
```

Dependency flow is one-directional. The coupler is pure grid/map machinery and
knows nothing about yelmox. `step_*` and `ice_domain` know the physics modules
but stay flavor-agnostic, so other yelmox variants can reuse them.

## Coupler (in fesm-utils, `coords` library)

Grids are identified by **string name** in every `remap` call. A name resolves to
its grid definition **from disk by default** — `grid_<name>.txt` in the map
folder, read via `grid_cdo_read_desc` (extended to parse both the fesm-utils `#`
header and cdo-native CF projection keys, so the existing `maps/grid_*.txt` are
used directly). An in-memory registry (`coupler_add_grid`) is an optional
override that wins when a name is registered. So nothing is hardcoded (no
`grid_select`): the available grids are whatever `grid_*.txt` files live in
`maps/`.

Fixed-capacity arrays + counters (not growable allocatables): appending to an
allocatable array of `map_class` would deep-copy every existing (large) map on
each add. Fixed shells cost nothing until `map_init` fills a map's allocatable
weight components.

```fortran
module coupler

    use coords, only : grid_class, map_class, map_init, map_field, grid_cdo_read_desc
    ! remap mirrors map_field's kind coverage (dp/sp/int); no wp<->dp conversion

    implicit none
    private

    integer, parameter :: GRID_MAX = 16, MAP_MAX = 64

    type grid_entry
        character(len=256) :: name
        type(grid_class)   :: grid          ! coords grid, what map_init consumes
    end type

    type coupler_class
        type(grid_entry) :: grids(GRID_MAX);  integer :: ngrids = 0
        type(map_class)  :: maps(MAP_MAX);     integer :: nmaps  = 0
        character(len=512) :: map_fldr = "maps"   ! disk cache dir for map_init
    end type

    interface remap
        ! mirror map_field's kind coverage; each wrapper allocates var_dst and
        ! calls the generic map_field (which dispatches by kind).
        module procedure remap_2d_dp, remap_2d_sp, remap_2d_int
        module procedure remap_3d_dp, remap_3d_sp, remap_3d_int
    end interface

    public :: coupler_class, coupler_init, coupler_add_grid, coupler_prime, remap

contains

    subroutine coupler_init(cpl, map_fldr)
        ! reset counters, set map_fldr
    end subroutine

    subroutine coupler_add_grid(cpl, name, grid)
        ! cpl%ngrids += 1; cpl%grids(cpl%ngrids) = grid_entry(name, grid)
        ! error on duplicate name or capacity overflow
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: name
        type(grid_class),    intent(in)    :: grid
    end subroutine

    subroutine coupler_prime(cpl, src, dst, method)
        ! eager warm-up: force-build a map up front (fail fast, cost up front).
        ! integer :: im; im = get_map(cpl, src, dst, method)   ! discard index
    end subroutine

    ! --- private: find-or-build the directed, method-typed map ---
    function get_map(cpl, src, dst, method) result(im)
        ! 1. linear-search cpl%maps(1:nmaps) for (name1==src, name2==dst, method)
        ! 2. on miss: resolve src/dst grids (in-memory registry, else read
        !    grid_<name>.txt from cpl%map_fldr); nmaps += 1;
        !    call map_init(cpl%maps(nmaps), grid_src, grid_dst, method=method,
        !                  fldr=cpl%map_fldr, load=.true.)   ! hits disk cache
        !    return nmaps
        integer :: im
    end function

    ! One wrapper per (kind x rank); dp shown, sp/int identical bar the type.
    subroutine remap_2d_dp(cpl, var_src, src, var_dst, dst, method, stat)
        type(coupler_class),   intent(inout) :: cpl
        real(dp),              intent(in)    :: var_src(:,:)
        character(len=*),      intent(in)    :: src, dst
        real(dp), allocatable, intent(inout) :: var_dst(:,:)     ! sized here
        character(len=*), optional, intent(in) :: method         ! default "con"
        character(len=*), optional, intent(in) :: stat           ! e.g. "mean"
        ! im = get_map(cpl, src, dst, method_or_default)
        ! allocate/re-shape var_dst to (nx_dst,ny_dst) if needed
        ! call map_field(cpl%maps(im), name, var_src, var_dst, stat=stat, method=method)
    end subroutine

    subroutine remap_3d_dp(cpl, var_src, src, var_dst, dst, method, stat)
        real(dp),              intent(in)    :: var_src(:,:,:)
        real(dp), allocatable, intent(inout) :: var_dst(:,:,:)   ! (nx,ny,nz_src)
        ! reshape var_dst to (nx_dst,ny_dst,size(var_src,3));
        ! loop trailing dim, apply the same 2D map per level
    end subroutine

end module coupler
```

`method` defaults to `"con"`; `"bilin"` is passed explicitly (CDO-consistent
spelling). `stat` is optional and forwarded to `map_field`.

## Domain + coupling steps (in `yelmox/libs`)

The whole state of one region is bundled into `ice_domain`. This makes the two
awkward requirements fall out naturally:

- **Coupling strategies** (which components are on/off): each component update is
  a reusable `step_*` primitive guarded internally by its `ctl` flag; the
  sequence is thin composition in `yelmox_step`.
- **Bipolar** (north + south): two explicit `type(ice_domain) :: dom_north,
  dom_south`. Each domain carries its own coupler/grids, so north and south never
  collide. Buffers are step-local allocatables (see below), inherently reentrant
  across domains. The single-domain driver runs the whole `yelmox_step`; the
  bipolar driver breaks it into its `step_*` primitives and interleaves them with
  the shared Ocean Box Model (see Drivers). Cross-hemisphere state — the
  barystatic sea level and the OBM — is driver-owned and shared, not per-domain.

```fortran
module yelmox_domain

    use coupler
    use yelmo_defs,    only : yelmo_class
    use marine_shelf,  only : marshelf_class, marshelf_update, marshelf_update_shelf
    use fastisostasy,  only : isos_class, isos_update
    use snapclim,      only : snapclim_class, snapclim_update
    use smbpal,        only : smbpal_class, smbpal_update_monthly
    use barysealevel,  only : bsl_class, bsl_update
    use htopo,         only : htopo_class, htopo_init

    type domain_ctl
        logical :: with_ice_sheet, with_isostasy, with_marine_shelf, with_climate
        character(len=256) :: domain       ! e.g. "Antarctica"
        character(len=256) :: grid_name    ! highest-res reference grid = htopo,
                                           !   e.g. "ANT-16KM" (ctl's top level)
        character(len=256) :: grid_yelmo   ! Yelmo grid, e.g. "ANT-32KM"
        character(len=256) :: grid_mshlf   ! marine-shelf grid (may equal grid_name)
        ! ... time control, output config
    end type

    type ice_domain
        type(yelmo_class)    :: yelmo
        type(marshelf_class) :: mshlf
        type(isos_class)     :: isos
        type(snapclim_class) :: snp
        type(smbpal_class)   :: smb
        ! NB: bsl (barystatic sea level) is NOT here -- it is a shared,
        !     driver-owned bsl_class, passed into the isostasy steps.
        type(htopo_class)    :: topo       ! hi-res geometry reference hub
        type(coupler_class)  :: cpl        ! this region's grids + maps
        type(domain_ctl)     :: ctl
    end type

contains

    subroutine domain_init(dom, path_par, ...)
        ! 1. init sub-models (Yelmo etc. on their own grids)
        ! 2. htopo_init(dom%topo, path_par, "htopo")   ! hi-res reference hub
        ! 3. coupler_init(dom%cpl)                      ! grids resolve from maps/*.txt
        ! 4. prime known maps (fail fast, cost up front):
        !      coupler_prime(dom%cpl, ctl%grid_yelmo, ctl%grid_name, "bilin")
        !      coupler_prime(dom%cpl, ctl%grid_name, ctl%grid_yelmo, "con")
        ! No coupler_add_grid needed: names resolve from grid_<name>.txt on disk.
    end subroutine

    subroutine yelmox_step(dom, time, bsl)  ! bsl: shared, driver-owned
        call step_isostasy(dom, time, bsl)  ! guarded by ctl%with_isostasy
        call step_icesheet(dom, time)       ! guarded by ctl%with_ice_sheet
        call step_climate(dom, time)        ! guarded by ctl%with_climate
        call step_marine_shelf(dom, time)   ! guarded by ctl%with_marine_shelf
    end subroutine

    subroutine step_marine_shelf(dom, time)
        real(wp), allocatable :: H_ice(:,:), z_bed(:,:), z_sl(:,:), f_grnd(:,:)
        real(wp), allocatable :: bmb(:,:), T_shlf(:,:)
        if (.not. dom%ctl%with_marine_shelf) return

        ! remap Yelmo -> mshlf grid (bilin)
        call remap(dom%cpl, dom%yelmo%tpo%now%H_ice, dom%ctl%grid_yelmo, &
                   H_ice, dom%ctl%grid_mshlf, method="bilin")
        call remap(dom%cpl, dom%yelmo%bnd%z_bed, dom%ctl%grid_yelmo, &
                   z_bed, dom%ctl%grid_mshlf, method="bilin")
        ! ... remaining inputs; unit-convert in place (newfield = f(buf))

        call marshelf_update_shelf(dom%mshlf, H_ice, z_bed, f_grnd, ..., dx=dx_mshlf)
        call marshelf_update(dom%mshlf, H_ice, z_bed, f_grnd, ..., dx=dx_mshlf)

        ! aggregate outputs back to Yelmo (con is the default)
        call remap(dom%cpl, dom%mshlf%now%bmb_shlf, dom%ctl%grid_mshlf, &
                   bmb, dom%ctl%grid_yelmo, stat="mean")
        call remap(dom%cpl, dom%mshlf%now%T_shlf, dom%ctl%grid_mshlf, &
                   T_shlf, dom%ctl%grid_yelmo, stat="mean")
        dom%yelmo%bnd%bmb_shlf = bmb
        dom%yelmo%bnd%T_shlf   = T_shlf
    end subroutine

    ! step_isostasy, step_climate, step_smb, step_icesheet — same shape
end module yelmox_domain
```

### Hi-res topography hub (htopo)

`htopo` sits *above* every physics module (including Yelmo): its grid is the
finest resolution in the setup, and it is the reference geometry the coupler
remaps *from*. On the topo grid it holds static masks `regions`/`basins` (loaded
once) and geometry `z_bed`/`H_ice`/`z_srf` (initial reference, later refreshed
each step from Yelmo/isostasy). It is configured by its own `[htopo]` namelist
group, whose `domain`/`grid_name` name the highest-res level and drive the
`{domain}/{grid_name}` path templating (`ctl%grid_name` mirrors `[htopo]
grid_name`):

```
&htopo
    domain       = "Antarctica"
    grid_name    = "ANT-16KM"
    topo_path    = "ice_data/{domain}/{grid_name}/{grid_name}_TOPO-BedMachine.nc"
    name_z_bed   = "z_bed"   name_H_ice = "H_ice"   name_z_srf = "z_srf"
    basins_path  = "ice_data/{domain}/{grid_name}/{grid_name}_BASINS-nasa.nc"
    name_basins  = "basin"
    regions_path = "ice_data/{domain}/{grid_name}/{grid_name}_REGIONS.nc"
    name_regions = "mask"
/
```

`htopo_init` resolves the grid from `grid_<name>.txt` (the disk grid table) and
reads the fields onto it — validated by `tests/test_htopo.f90` against the real
ANT-16KM data.

### Buffers

Fields that cross a grid boundary land in **step-local allocatables** (e.g.
`H_ice` above), not a coupler-owned pool. `remap` takes the destination as
`intent(inout), allocatable` and sizes it itself (allocate-on-demand, reshape if
wrong), so callers never hand-compute `nx,ny`. Unit conversions are done in place
on the buffer (`newfield = f(buf)`); they are caller concerns, kept out of the
coupler. 3D (e.g. monthly) fields use the `remap_3d` overload. Step-local
allocatables are reentrant, which is what bipolar needs; per-timestep
reallocation cost is negligible against the physics.

## Drivers

Three thin programs share `yelmox_domain` (all per-domain physics/coupling lives
there, so the drivers only differ in config parsing + the loop over domains):

- **`yelmox`** (single domain) — argument is one domain nml; one `ice_domain`,
  output to the run dir.
- **`yelmox_bipolar`** (bipolar, in `yelmox_bipolar/`) — argument is a single parameter
  file holding both hemispheres, in the original `yelmox_bipolar` convention: each
  domain's instance groups carry a hemisphere suffix (`yelmo_south`,
  `coupling_north`, `snap_south`, …), while `[ctrl]`, `[barysealevel]`, the
  `[nautilus]`/`[stommel]` OBM block and the yelmo physics groups (`ydyn`,
  `ytopo`, …) are shared. `[ctrl] active_north`/`active_south` select the
  hemispheres (hardcoded `_north`/`_south` suffixes threaded into every `group=`
  via `domain_init(..., group_suffix=…)`). Distinct group names also let
  `runme -p group.name=val` target one hemisphere unambiguously. Each domain
  writes to a subfolder named after its domain; invoke with
  `runme -e bipolar -n yelmox_bipolar/yelmox_bipolar_Bipolar.nml`.

  A bipolar run is never more than north + south, so the two domains are explicit
  variables, not an array — the ocean coupling is asymmetric (north ↔
  `obm%fn/thetan/tn`, south ↔ `obm%fs/thetas/ts`). The driver owns the shared
  `bsl` and `obm`, breaks `yelmox_step` into its `step_*` primitives, and
  interleaves the OBM exactly as in `yelmox_bipolar` (below). The ocean exchanges
  live in `yelmox_osm/obm_coupling.f90`; they and the OBM default off (`[ctrl]
  active_obm=False`), so the config runs as two independent domains until enabled.

```fortran
program yelmox_bipolar
    use yelmox_domain
    use obm, only : obm_update
    use obm_coupling
    type(ice_domain) :: dom_north, dom_south
    type(bsl_class)  :: bsl        ! shared, driver-owned
    type(obm_class)  :: obm        ! shared, driver-owned

    ! read one file -> shared timeline + active_north/active_south
    if (active_north) call domain_init(dom_north, path_par, ..., group_suffix="_north")
    if (active_south) call domain_init(dom_south, path_par, ..., group_suffix="_south")

    do while (time < time_end)
        call bsl_update(bsl, ...)                          ! once, shared
        call advance_isostasy(dom_north); call advance_isostasy(dom_south)
        if (active_obm) call obm_update(obm, dtt, obm_name)
        call advance_dynamics(dom_north); call advance_dynamics(dom_south)
        ! ocean coupling (per hemisphere): atm2obm, ism2obm (fwf), hyster, obm2ism
        call step_marine_shelf(dom_north, ...); call step_marine_shelf(dom_south, ...)
        ! per-domain output/restart (subfolder each) + shared bsl/obm restart
        time = time + dtt
    end do
end program
```

- **`yelmox_esm`** (single domain, in `yelmox_esm/`) — ESM climatic forcing in
  place of snapclim. Reuses `domain_init` (with `init_climate=.false.`, so
  snapclim is skipped) plus the shared `step_optimize/step_isostasy/step_icesheet/
  refresh_htopo` primitives and the restart bundle, but the driver owns an
  `esm_forcing_class` and calls its own `step_climate_esm` / `step_marine_shelf_esm`
  (contained in the program) instead of the snapclim-based steps. esm is a
  first-class multigrid component: it runs entirely on its own grid (*esm grid* =
  `grid_clim`, exactly like snapclim), and each output is remapped to the consumer
  module's grid at coupling time — atmosphere to `grid_smb` (smbpal), the
  depth-interpolated ocean forcing to `grid_mshlf` (marine_shelf). This works
  because `marshelf_interp_shelf` reads only `mshlf%par` (grid-agnostic), so the ESM
  ocean interpolation runs on `grid_clim` and the resulting `T_shlf`/`S_shlf` are
  remapped to `grid_mshlf` before `marshelf_update`; `esm.f90` is untouched.
  Geometry comes from the hub, remapped to whichever grid a step needs; SMB / ocean
  BCs aggregate back to Yelmo. With `grid_clim == grid_smb == grid_mshlf ==
  grid_name == grid_yelmo` every remap is an identity copy, reproducing
  `yelmox_esm.f90`; set `grid_clim` to a coarse ESM grid and it genuinely fans out.
  Config splits ESM-specific control ([esm] + the run_step group
  [spinup]/[transient]: `time_ref/hist/proj/esm_ref`, `use_*`, CMIP output) from
  the shared mg groups ([coupling]/[output]/[htopo]). Output (incl. the
  CMIP-formatted files) is kept identical to `yelmox_esm.f90` via the
  `yelmox_esm_output` module (in the same folder). Invoke with
  `runme -e esm -n yelmox_esm/yelmox_esm_Antarctica.nml`.

```fortran
program yelmox_esm
    use yelmox_domain
    use esm
    use yelmox_esm_output
    type(ice_domain)        :: dom
    type(bsl_class)         :: bsl    ! shared, driver-owned
    type(esm_forcing_class) :: esm    ! driver-owned climate (replaces snapclim)

    call domain_init(dom, path_par, ts%time, init_climate=.false.)   ! skip snapclim
    call esm_forcing_init(esm, ..., grid_name=dom%ctl%grid_clim)      ! on esm's own grid

    do while (.not. ts%is_finished)
        call bsl_update(bsl, ...)                       ! once, shared
        call step_optimize(dom, ts); call step_isostasy(dom, ts, bsl)
        call step_icesheet(dom, ts); call refresh_htopo(dom)
        call step_climate_esm(dom, esm, ec, ts)         ! esm + smbpal (contained)
        call step_marine_shelf_esm(dom, esm, ec, ts)    ! esm ocean BCs (contained)
        ! output (yelmo2D / yelmo1D_esm / CMIP) + restart bundle + shared bsl
    end do
end program
```

## Integration gaps to resolve during the build

1. **Grid definitions come from disk.** `grid_<name>.txt` in `maps/` (parsed by
   `grid_cdo_read_desc`, incl. cdo-native CF keys) — no `grid_select`, no
   `ygrid_class`→`grid_class` conversion. Assumes `grid_<yelmo>.txt` matches
   Yelmo's runtime grid (true for standard full-domain runs); the in-memory
   registry override is the escape hatch otherwise. *(Done.)*
2. **Precision.** No external conversion needed — `map_field` already has
   `dp`/`sp`/`int` variants. `remap` provides matching `dp`/`sp`/`int` × 2d/3d
   wrappers under one generic interface, each calling the generic `map_field`.
   *(Done.)*
3. **Grid names** are sourced from `domain_ctl` (and `[htopo]`) as the source of
   truth for `remap` keys.
4. **Conservative area basis.** `"con"` weights on projected cell area — sanity
   check mass conservation of `z_bed`/`bmb` aggregation on a real grid pair before
   trusting it as a BC.
5. **FastIsostasy hi-res output** — deferred. Lean toward the coupler upscaling
   the 16KM isostasy output rather than making the solver grid-aware. Not in the
   first skeleton.

## Commit order

1. `coupler.f90` in fesm-utils + `test_coupler` (remap 16↔2KM both directions,
   `con` conservation, cache/prime, disk-driven resolution). *(Done.)*
2. `grid_cdo_read_desc` cdo-native CF-key parsing + `test_grid_cf_read`;
   coupler disk resolution of `grid_<name>.txt`. *(Done.)*
3. `yelmox_domain.f90` skeleton: `ice_domain`, `domain_ctl`, empty `step_*`. *(Done.)*
4. `htopo.f90` hi-res reference hub + `test_htopo` (load ANT-16KM fields). *(Done.)*
5. `domain_init` (sub-model init on their grids + htopo + map priming). *(Done.)*
6. Fill `step_*` one module at a time (marine_shelf first); diff vs `yelmox.f90`. *(Done: isostasy, ice sheet, climate, smb, marine_shelf, optimization.)*
7. `yelmox.f90` driver; validate single-domain parity with `yelmox.f90`. *(Done: identity-grid parity tracks the reference; the residual is a small deterministic startup offset, independent of isostasy — bit parity not required for the move.)*
8. Enable `nd=2` bipolar. *(Done: driver takes one par file per domain; AIS+GRL runs on a shared timeline to per-domain subfolders.)*

Also ported from `yelmox.f90`: `equil_method="opt"` (basal-friction + thermal-forcing
optimization, `step_optimize`), `smb_method="smb_simple"` + `calc_glacial_smb`, and
domain-specific cold-start setup (Antarctica equilibration, Greenland masks/marine-ice/
NEGIS, Laurentide/North LGM initialization).
