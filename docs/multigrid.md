# Multigrid coupling

Architecture of the two YelmoX programs, `yelmox` and `yelmox_bipolar`. Each
component of a domain (Yelmo, isostasy, climate, smb, marine shelf) runs on its
own grid, and fields are remapped between grids at the moment of coupling. The
domain is defined around a hi-resolution geometry hub (htopo), the finest grid of
the setup, from which every component, Yelmo included, is populated.

Example for Antarctica with Yelmo on `ANT-16KM` and a hi-res hub on `ANT-2KM`:

- `H_ice` (Yelmo `ANT-16KM`) → bilin → `ANT-2KM` for the marine-shelf calc.
- `z_bed`/`z_ss` from FastIsostasy → conservatively aggregated to Yelmo's grid,
  set as BCs.
- marine_shelf computed at `ANT-2KM`; `bmb_shlf`/`T_shlf` → conservatively
  aggregated → `ANT-16KM` → passed to Yelmo as forcing.

The bookkeeping of the remapping steps lives in a coupler (fesm-utils) and in
the domain's coupling steps (yelmox), not in the component modules.

## Building blocks

- **Component modules are grid-agnostic.** `marshelf_update`,
  `isos_update`, `snapclim_update`, `smbpal_update` all take/return plain
  `(nx,ny)` arrays; the grid is the caller's choice. Moving a module onto
  another grid is a caller-side change, not a module change.
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
fesm-utils/src/coupler.f90        coupler_class + remap (grids, map cache)
        │
        ▼
yelmox/libs/
  kryos.f90             kryos_domain, domain_ctl, domain_init, remap
  kryos_coupling.f90    step_* and couple_* primitives
  kryos_startup.f90     cold start, restart bundles (domain_startup, ...)
  kryos_output.f90      per-module 2D/1D output, CMIP output
  kryos_regions.f90     named regions, marine-ice start, NEGIS, glacial smb
  kryos_forcing.f90     driver-owned transient forcing (tsgen)
  htopo.f90             the hi-res geometry hub
  yelmox_climate.f90    climate backends behind one interface (climate_out)
        │
        ▼
yelmox/yelmox.f90                 single domain
yelmox_bipolar/yelmox_bipolar.f90 north + south, shared sea level + OBM
```

Dependency flow is one-directional. The coupler is pure grid/map machinery and
knows nothing about yelmox. The Kryos modules know the physics modules but no
driver, so both programs share them.

## Coupler (fesm-utils)

```fortran
call coupler_init(cpl, map_fldr, gen)               ! map folder, default map generator
call coupler_add_grid(cpl, name, grid)              ! optional in-memory grid
call coupler_prime(cpl, src, dst, method, gen)      ! build a map up front
call remap(cpl, var_src, src, var_dst, dst, method, stat)
```

Grids are identified by **string name** in every `remap` call. A name resolves to
its grid definition **from disk by default** — `grid_<name>.txt` in the map
folder (`maps/`), read via `grid_cdo_read_desc` (both the fesm-utils `#` header
and cdo-native CF projection keys). `coupler_add_grid` registers an in-memory
grid that wins over the disk definition. Nothing is hardcoded: the available
grids are the `grid_*.txt` files in `maps/`.

`remap` finds the map for `(src, dst, method)` in the cache, or builds it (from
its disk cache when present) and stores it. `method` is `"con"` (conservative,
the default) or a distance kernel (`"nn"`, `"shepard"`, `"quadrant"`,
`"bilin"`); `stat` (`mean`/`count`/`stdev`) is forwarded per call. `gen`
(`"coords"`, the default, or `"cdo"`) only selects how the weights are built.
`remap` covers `dp`/`sp`/`int` × 2D/3D, and sizes the destination itself. The
maps live in fixed-capacity storage, so adding a map never deep-copies the
(large) maps already built.

The domain wraps it as `remap(dom, var_src, src, var_dst, dst, method)`
(`kryos`), which copies when `src == dst`, so a component on Yelmo's grid costs
nothing.

## The domain (`kryos_domain`)

The whole state of one region is bundled into `kryos_domain`:

```fortran
type kryos_domain
    type(yelmo_class)          :: yelmo
    type(marshelf_class)       :: mshlf
    type(isos_class)           :: isos
    type(yelmox_climate_class) :: cl     ! climate backend ([comps] climate)
    type(climate_out_class)    :: clim   ! backend-agnostic climate output (now/ref)
    type(smbpal_class)         :: smb
    type(smb_simple_class)     :: smbs   ! surface_method = "smb_simple"
    type(sediments_class)      :: sed
    type(geothermal_class)     :: gthrm
    type(phys_const_class)     :: cnst   ! physical constants, shared by every component
    type(htopo_class)          :: topo   ! hi-res geometry hub
    type(coupler_class)        :: cpl    ! this region's grids + map cache
    type(ice_opt_params)       :: opt    ! spin-up optimization
    type(negis_params)         :: ngs
    type(glacial_smb_params)   :: gsmb
    type(domain_ctl)           :: ctl    ! [domain], [comps], [sim], [output]
end type
```

The barystatic sea level (`bsl`) and the transient forcing (`tsf`) are not part
of the domain: they are driver-owned, so a multi-domain driver shares them (or
not). Each `step_*` primitive advances one component on its own grid and is a
no-op when that component is inactive (`with_*`); each `couple_*_to_yelmo` lands
one component's output on Yelmo's grid. Two domains (bipolar) carry their own
coupler and grids, so they never collide.

The climate reaches the rest of the domain only through `climate_out_class`
(`now` and `ref`): the atmosphere, the ocean as depth profiles or at the shelf
base, and, when the backend has them, the surface mass balance and subglacial
discharge. Nothing downstream depends on which backend produced them.

### Domain definition (`[domain]`) and the hi-res hub (htopo)

One `[domain]` group (`[domain_north]`/`[domain_south]` in bipolar) defines the
domain: its name, its physical constants, the grid of every component, and the
hub's topography and code masks. Yelmo takes the domain name and its grid from it
(`[yelmo]` no longer sets `domain`/`grid_name`). A blank component grid takes its
default:

```
&domain
    name         = "Antarctica"
    phys_const   = "Earth"      ! physical constants: group of input/yelmo_phys_const.nml
    grid_hub     = "ANT-16KM"   ! hi-res geometry hub
    grid_ice     = "ANT-32KM"   ! Yelmo                                [grid_hub]
    grid_isos    = ""           ! isostasy                             [grid_ice]
    grid_clim    = ""           ! reference climate + transient forcing [grid_ice]
    grid_surface = ""           ! surface mass balance                 [grid_clim]
    grid_shelf   = ""           ! marine shelf                         [grid_hub]
    topo_path    = "ice_data/{domain}/{grid_name}/{grid_name}_TOPO-BedMachine.nc"
    topo_names   = "z_bed" "H_ice" "z_srf" "z_bed_sd"   ! z_bed_sd: "" = none (0)
    regions_path = "ice_data/{domain}/{grid_name}/{grid_name}_REGIONS.nc"   ! "" = none (1)
    regions_var  = "mask"
    basins_path  = "ice_data/{domain}/{grid_name}/{grid_name}_BASINS-nasa.nc"
    basins_var   = "basin"
    sectors_path = "ice_data/{domain}/{grid_name}/{grid_name}_BASINS-nasa.nc"   ! "" = none (1)
    sectors_var  = "mask_regions"
    ice_codes_mode = "exclude"      ! where ice is allowed: all | include | exclude (ice_codes)
    ice_codes    = 2.0              ! codes of regions ("" = none)
    region_names = "APIS" "WAIS" "EAIS"   ! named regions for 1D output ("" = none)
    region_mask  = "sectors"        ! regions | basins | sectors
    region_codes = 3.0 1.0 2.0
    relax_codes_mode = "none"       ! where ice relaxes to the reference: none | all | include | exclude
    relax_codes  = ""               ! codes of regions ("" = none)
    relax_tau    = 0.0              ! [yr] relaxation timescale there
/
```

The domain loads its physical constants once (`phys_const_load`, the group
`phys_const` of `input/yelmo_phys_const.nml`) and hands the same record to every
component: the hub, Yelmo (`yelmo_init` `cnst`), isostasy, the marine shelf and
`smb_simple`. Yelmo's own `yelmo.phys_const` then only selects its calendar year
(`sec_year`).

`htopo` holds the hub. It sits *above* every physics module (including Yelmo):
its grid (`grid_hub`) is the finest resolution in the setup, and it is the
reference geometry the coupler remaps *from*. On the hub grid it holds static
fields (the code masks `regions`/`basins`/`sectors` and the bed roughness
`z_bed_sd`, loaded once), the reference geometry `z_bed_ref`/`H_ice_ref`/`z_srf_ref`
(loaded once) and the current geometry `z_bed`/`H_ice`/`z_srf`/`f_grnd`/`z_sl`
(refreshed each step, see below). `{domain}/{grid_name}` in the paths
resolve to `name`/`grid_hub`. `htopo_init` resolves the grid from
`grid_<name>.txt` (the disk grid table) and reads the fields onto it —
validated by `tests/test_htopo.f90` against the real ANT-16KM data. A blank
mask path (e.g. paleo domains without mask files) leaves the mask at `1.0`
(single region/basin/sector); a blank `z_bed_sd` name leaves it at `0`. Gaps in the
topography file (missing values, e.g. outside the coverage of the source
dataset) are filled: no ice, the bed from the nearest valid cell, the surface
from the bed and the ice thickness (sea level 0), and no bed roughness.

Yelmo is populated from the domain like the other components. Its grid comes
from `maps/grid_<grid_ice>.txt` (`yelmo_init_grid`, `grid_def="none"`), and the
hub's `z_bed`/`H_ice`/`z_srf`/`z_bed_sd`, remapped conservatively to `grid_ice`,
are both its initial topography and its present-day reference (`yelmo_init`
`topo_init`/`topo_pd`). Yelmo then processes them as it would its own files
(`[yelmo_init_topo]` keeps `init_topo_state`, `z_bed_f_sd`, smoothing; its
`grad_lim_zb` applies). Where `grid_ice = grid_hub` the remap is a copy.

The hub follows the models each step (`couple_yelmo_to_htopo`, after Yelmo). On
Yelmo's grid it mirrors Yelmo, fractional grounding included. On a finer hub it keeps
its hi-res reference and adds Yelmo's anomalies, refined bilinearly: the bed
displacement `z_bed - z_bed_ref`, and the change in ice thickness from the hub
reference as Yelmo received it (remapped conservatively), clipped at zero
thickness. `htopo_update` then recomputes on the hub the grounding (0 or 1, from
flotation) and the surface elevation, with Yelmo's densities; sea level is
refined bilinearly. The hi-res bed and ice therefore reach the marine shelf and
the climate's surface elevation, instead of a refined copy of Yelmo's fields.

The code masks reach Yelmo the same way (nearest neighbour): `regions` and
`basins` are Yelmo's (`yelmo_init` `regions`/`basins`), and every component uses
this one set. Where ice is allowed follows from `ice_codes_mode` and
`ice_codes` (codes of `regions`; `yelmo_init` `mask_ice`); Yelmo's
`mask_border` (`[yelmo]`, default `"auto"`) then sets the domain border.
Where the ice relaxes to the reference follows from `relax_codes_mode` and
`relax_codes` (codes of `regions`): Yelmo's `tau_relax` is `relax_tau` there and
-1 (free) elsewhere, used with `ytopo.topo_rel = -1`. The
named regions (`region_names`, one code each of `region_mask`) get their own 1D
output, `yelmo_ts_<name>.nc`.

### Buffers

Fields that cross a grid boundary land in **step-local allocatables**, not a
coupler-owned pool. `remap` takes the destination as
`intent(inout), allocatable` and sizes it itself (allocate-on-demand, reshape if
wrong), so callers never hand-compute `nx,ny`. Unit conversions are done in place
on the buffer (`newfield = f(buf)`); they are caller concerns, kept out of the
coupler. 3D (e.g. monthly) fields use the `remap_3d` overload. Step-local
allocatables are reentrant, which is what bipolar needs; per-timestep
reallocation cost is negligible against the physics.

## Drivers

Both programs use the same driver plumbing:

- **`tstep_init(ts, path_par, group, dtt [, time_ref, cal])`** (fesm-utils
  `timestepping`) — reads the run's timeline group (`[ctrl] run_step`: `"ctrl"`
  for `[ctrl]` itself, or a phase group such as `[spinup]`/`[transient]`) and
  initializes the driver-owned timestepper; `[ctrl] calendar` / `calendar_ref`
  give `cal` / `time_ref` (calendar years, e.g. ESM runs). `domain_init` takes
  the same group name (`timeline_group`) and reads the values the domain needs
  (`tstep_method`, `dtt`) itself.
- **`domain_startup(dom, ts, bsl, restart [, restore_bsl, tsf])`**
  (`kryos_startup`) — cold start (`restart = "None"`, `domain_init_state`) or
  restore of the domain bundle `restart` + hub rebuild. The driver reads
  `[ctrl] restart`. `yelmox` restores the shared bsl from the same bundle;
  `yelmox_bipolar` restores it once via **`bsl_startup(bsl, fldr)`**, passes
  each domain its subfolder of the bundle and `restore_bsl=.false.`.
- **`domain_init_ice(dom, ts)`** — the cold-start ice state after
  `yelmo_init_state` (`[sim]` `init_kill_shelves`, `init_marine_H`,
  `init_method`, `init_time_thrm`; see [yelmox](flavor-yelmox.md#configuration)).
- **`run_restart_write(dom, bsl, time [, tsf])`** — the single-domain restart
  bundle (domain sub-models + `bsl_restart.nc` + the tsforcing state, one
  auto-named folder).
- **`cadence_due(time, dt)`** (`kryos`) — the cadence predicate for `dt_clim` and
  `dt_cmip`; `dt <= 0` disables a cadence. Output and restart bundles follow
  fesm-utils `timeout` schedules (`[tm_1D]`, `[tm_2D]`, `[tm_2Dsm]`, `[tm_rst]`);
  a restart bundle is always written at `time_end`.

Every driver loop has the same shape: output and restarts are written at the
top of the loop for the current time (`time_init` on the first pass), then the
loop exits if the run is finished (`time_end`, or a tripped kill switch), else
`tstep_update` advances the time and the domain is stepped with the coupling
sequence written out inline. Each output call appears once, and the final state
+ restart bundle are written on the last pass (`timeout_check(...) .or.
ts%is_finished`).

- **`yelmox`** — argument is one parameter file; one `kryos_domain`, output to
  the run dir. The climate backend is chosen at runtime (`[comps] climate`),
  so ESM and REMBO runs use this program too. See [yelmox](flavor-yelmox.md).
- **`yelmox_bipolar`** — argument is one parameter file holding both
  hemispheres. Each domain's groups carry a hemisphere suffix (`yelmo_south`,
  `domain_north`, `comps_north`, `snap_south`, …), threaded into every group
  via `domain_init(..., group_suffix=…)`; `[ctrl]`, `[barysealevel]`, the OBM
  groups and the Yelmo physics groups (`ydyn`, `ytopo`, …) are shared. Distinct
  group names also let `runme -p group.name=val` target one hemisphere. The two
  domains are explicit variables (`dom_north`, `dom_south`), not an array: the
  ocean coupling is asymmetric. The driver owns the shared `bsl` and Ocean Box
  Model and interleaves them with the per-domain steps; the ocean coupling lives
  in `yelmox_bipolar/obm_coupling.f90`. Each domain writes to a subfolder named
  after it. See [yelmox_bipolar](flavor-bipolar.md).

## Open issues

1. **Conservative area basis.** `"con"` weights use projected cell area — sanity
   check mass conservation of `z_bed`/`bmb` aggregation on a real grid pair before
   trusting it as a BC.
2. **FastIsostasy hi-res output** — deferred. Lean toward the coupler refining
   the isostasy output rather than making the solver grid-aware.
3. **One grid for reference climate and transient forcing.** `grid_clim` sets
   the grid of both the reference climatology (often from a high-resolution
   regional model) and the transient forcing (often from a coarser climate
   model). A coarse `grid_clim` matches the forcing but loses the detail of the
   high-res reference. Giving the two their own grids is a candidate for future
   work; until then, set `grid_clim` to the highest-resolution climate input.
