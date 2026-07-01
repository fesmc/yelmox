# Multigrid coupling (`yelmox_mg`)

Design doc for a multigrid rewrite of the yelmox driver. Status: **design, pre-implementation.**

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
yelmox/yelmox_mg.f90                        thin driver: init → time loop → I/O
```

Dependency flow is one-directional. The coupler is pure grid/map machinery and
knows nothing about yelmox. `step_*` and `ice_domain` know the physics modules
but stay flavor-agnostic, so other yelmox variants can reuse them.

## Coupler (in fesm-utils, `coords` library)

Grids are identified by **string name** in every `remap` call — the name is the
key the coupler uses to look up the registered grid. Callers pass names sourced
from config/objects, never grid objects (the grids already live in the coupler).

Fixed-capacity arrays + counters (not growable allocatables): appending to an
allocatable array of `map_class` would deep-copy every existing (large) map on
each add. Fixed shells cost nothing until `map_init` fills a map's allocatable
weight components.

```fortran
module coupler

    use coords, only : grid_class, grid_init, map_class, map_init, map_field
    ! precision: map_field is dp; wrappers convert wp<->dp

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
        ! 2. on miss: find src/dst grids in registry by name; nmaps += 1;
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
- **Bipolar** (north + south): `type(ice_domain) :: dom(2)`, loop over them
  through the same `yelmox_step`. Each domain carries its own coupler/grids, so
  north and south never collide. Buffers are step-local allocatables (see below),
  which are inherently reentrant across domains.

```fortran
module yelmox_domain

    use coupler
    use yelmo_defs,    only : yelmo_class
    use marine_shelf,  only : marshelf_class, marshelf_update, marshelf_update_shelf
    use fastisostasy,  only : isos_class, isos_update
    use snapclim,      only : snapclim_class, snapclim_update
    use smbpal,        only : smbpal_class, smbpal_update_monthly
    use barysealevel,  only : bsl_class, bsl_update

    type domain_ctl
        logical :: with_ice_sheet, with_isostasy, with_marine_shelf, with_climate
        character(len=256) :: grid_yelmo   ! e.g. "ANT-16KM"  (source of truth
        character(len=256) :: grid_mshlf   ! e.g. "ANT-2KM"    for grid names)
        character(len=256) :: grid_topo    ! hi-res topo grid, loaded from file
        ! ... time control, output config
    end type

    type ice_domain
        type(yelmo_class)    :: yelmo
        type(marshelf_class) :: mshlf
        type(isos_class)     :: isos
        type(snapclim_class) :: snp
        type(smbpal_class)   :: smb
        type(bsl_class)      :: bsl
        type(coupler_class)  :: cpl        ! this region's grids + maps
        type(domain_ctl)     :: ctl
    end type

contains

    subroutine domain_init(dom, path_par, ...)
        ! 1. init sub-models
        ! 2. coupler_init(dom%cpl)
        ! 3. register grids:
        !      coupler_add_grid(dom%cpl, ctl%grid_yelmo, <coords grid from ygrid>)
        !      coupler_add_grid(dom%cpl, ctl%grid_topo,  <coords grid from file>)
        !      coupler_add_grid(dom%cpl, ctl%grid_mshlf, ...)
        ! 4. prime known maps (fail fast, cost up front):
        !      coupler_prime(dom%cpl, ctl%grid_yelmo, ctl%grid_mshlf, "bilin")
        !      coupler_prime(dom%cpl, ctl%grid_mshlf, ctl%grid_yelmo, "con")
    end subroutine

    subroutine yelmox_step(dom, time)
        call step_isostasy(dom, time)       ! guarded by ctl%with_isostasy
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

### Buffers

Fields that cross a grid boundary land in **step-local allocatables** (e.g.
`H_ice` above), not a coupler-owned pool. `remap` takes the destination as
`intent(inout), allocatable` and sizes it itself (allocate-on-demand, reshape if
wrong), so callers never hand-compute `nx,ny`. Unit conversions are done in place
on the buffer (`newfield = f(buf)`); they are caller concerns, kept out of the
coupler. 3D (e.g. monthly) fields use the `remap_3d` overload. Step-local
allocatables are reentrant, which is what bipolar needs; per-timestep
reallocation cost is negligible against the physics.

## Driver (`yelmox_mg.f90`)

```fortran
program yelmox_mg
    use yelmox_domain
    type(ice_domain), allocatable :: dom(:)      ! (1) single, (2) bipolar
    integer :: nd, k

    ! read config -> nd + per-domain path_par
    allocate(dom(nd))
    do k = 1, nd
        call domain_init(dom(k), path_par(k), ...)
    end do

    do while (time < time_end)
        do k = 1, nd
            call yelmox_step(dom(k), time)
        end do
        ! per-domain restart + output writing
        time = time + dt
    end do
end program
```

## Integration gaps to resolve during the build

1. **`ygrid_class` → coords `grid_class`.** `map_init` needs a coords grid; Yelmo
   exposes `ygrid_class`. Build the coords grid via `grid_init` from Yelmo's
   projection params + `xc/yc`. Verify projection metadata lines up against a real
   Yelmo grid before committing.
2. **Precision.** No external conversion needed — `map_field` already has
   `dp`/`sp`/`int` variants. `remap` provides matching `dp`/`sp`/`int` × 2d/3d
   wrappers under one generic interface, each calling the generic `map_field`.
3. **`marshelf_grid_class` gains a `name` field** (agreed) for consistency; grid
   names remain sourced from `domain_ctl` as the source of truth.
4. **Conservative area basis.** `"con"` weights on projected cell area — sanity
   check mass conservation of `z_bed`/`bmb` aggregation on a real grid pair before
   trusting it as a BC.
5. **FastIsostasy hi-res output** — deferred. Lean toward the coupler upscaling
   the 16KM isostasy output rather than making the solver grid-aware. Not in the
   first skeleton.

## Commit order (worktree)

1. `coupler.f90` in fesm-utils + a standalone unit test (remap a known field
   16↔2KM both directions; check `con` conservation).
2. `yelmox_domain.f90` skeleton: `ice_domain`, `domain_ctl`, empty `step_*`.
3. Fill `step_*` one module at a time; diff behavior against `yelmox.f90` on a
   single domain.
4. `yelmox_mg.f90` driver; validate single-domain parity with `yelmox.f90`.
5. Enable `nd=2` bipolar.
