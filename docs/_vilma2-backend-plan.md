# VILMA2 as an isostasy backend: plan

Status: design (2026-10-23). The leading `_` keeps quarto from rendering it.

States checked: yelmox dev a0be9e5 (v2.4 + stable pins), FastIsostasy (yelmox
submodule), fesm-utils dev d994f08 (v1.4), vilma2 60a85df, FesmData (local).

## Goal

Switch at runtime between FastIsostasy (default) and VILMA2 for bedrock and
sea level:

```
[isostasy]
    method = "fastisostasy"    ! "fastisostasy" | "vilma2"
```

VILMA2 is optional at build time (`vilma=1`), with a stub otherwise, as for
rembo. VILMA2 itself is still "under heavy development"; yelmox keeps
FastIsostasy as the default and treats VILMA2 as experimental until it is
validated in coupled runs.

## Decisions

- New module `libs/kryos_isostasy.f90`: the only place kryos touches an
  isostasy backend.
- Grids: no special handling. VILMA2 runs on its own grid like any other
  component; the coupler (fesm-utils `coords`) remaps. `grid_cdo` already
  reads `gaussian` and `lonlat` descriptions.
- VILMA2 is extended to return bed uplift and sea-surface change separately
  (not only rsl), and a barystatic sea level.
- Under `vilma2`, bsl comes from VILMA2. VILMA2 always sees a global load.
- Ice outside the domains: both a fixed present-day reference and a global
  ice history from file (`ice_ext`).
- Global reference topography and the 3D viscosity file are namelist
  options; the defaults are FesmData products (topography, Bagge2021_visc).
- `dzbdt_corr` has no VILMA2 analogue: ignored under `vilma2`, with a
  warning when it is non-zero.
- The full `&vilma` group is configurable from the yelmox par file. Its
  defaults are `input/vilma_defaults.nml`, a synchronized copy of VILMA2's
  (the chion/elsa convention).
- Build: `vilma=1` links VILMA2; default builds compile a stub.

## Current coupling (FastIsostasy)

Per domain, `dom%isos` is a `fastisostasy::isos_class` on `grid_isos`
(`libs/kryos.f90:166`). The driver owns one `bsl_class`, shared by all
domains (`yelmox/yelmox.f90:38`, `yelmox_bipolar/yelmox_bipolar.f90:64`).

| Call site | What |
|---|---|
| `kryos.f90:348` | `isos_init(dom%isos, ..., nx_i, ny_i, dx_isos, dy_isos, ...)` |
| `kryos_startup.f90:218` `domain_init_isostasy` | `isos_init_ref` + `isos_init_state` (cold start) |
| `kryos_startup.f90:~500` `domain_restart_read` | `isos_init_state` with `use_restart` |
| `kryos_startup.f90:421` | `isos_restart_write` |
| `kryos_coupling.f90:119` `step_isostasy` | remap `H_ice`, `dzbdt_corr` to `grid_isos`; `isos_update(..., bsl, dwdt_corr)` |
| `kryos_coupling.f90:186` `couple_isostasy_to_yelmo` | `z_bed = z_bed_ref + remap(w + we)`; `z_sl = bsl + remap(z_ss - bsl)` |
| `kryos_coupling.f90:229` `check_isostasy_reference` | reference bed vs Yelmo `z_bed_ref` on `grid_isos` |
| `kryos_output.f90:260,383` | `isos.nc` (z_bed, z_ss, w, we), `isos_ts.nc` (bsl, means) |
| drivers | `bsl_init`, `bsl_update`, `bsl_restart_write`, `bsl_startup` |

Only anomalies cross into Yelmo (`dz_bed`, `dz_ss`), so Yelmo keeps its own
high-resolution reference bed. That contract carries over unchanged.

## VILMA2 host API (as is)

`use vilma`, type `solid_earth` (`src/vilma_coupling.f90`):

- `vilma_par_load(se%par, filename, defaults_file, group)`: `&vilma` group
  overlaid on a complete defaults file (yelmo convention).
- `solid_earth_init(se, z_bed_eq, h_ice_eq, grid, h_ice_init)`: global
  reference bed and grounded ice on the host grid. If `grid` is the model's
  Gauss grid (`grid_is_gauss`), fields pass straight through; otherwise
  VILMA2 builds its own lon-lat → Gauss map.
- `solid_earth_spinup(se, h_ice_lgm)`: optional LGM-memory spin-up.
- `solid_earth_update(se, h_ice, dt_yr)`: advance by `dt_yr`; the load is
  ramped linearly over the interval. Returns `se%rsl`,
  `se%z_bed = z_bed_eq - rsl`, `se%bsl`.
- `vilma_restart_write/read`, `solid_earth_finalize`.
- Precision `wp = real64`; depends on fesm-utils `coords`, FFTW and SHTns.

Gaps for yelmox:

1. Only `rsl = N - u + Δφ` reaches the host. Yelmo needs bed uplift and
   sea-surface change separately (`z_bed`, `z_sl`). `u_lm`/`N_lm` exist
   inside the SLE (`vilma_sle.f90:55`) but are not exported.
2. `se%bsl` is an ice-volume-equivalent diagnostic, not the barystatic rise
   (see the comment on `update_bsl`, `vilma_coupling.f90:573`).
3. Only executables are built; there is no `libvilma.a` target.
4. The README asks for fesm-utils `coords-dev`; yelmox's fesm-utils dev
   already carries `coords` (coords-tmerc merged). Needs a build check
   against fesm-utils dev.
5. The 3D viscosity reader expects the isostasy_data file
   (`input/bagge2021.nc`: one structure, `lgvisc` on `radius`). FesmData's
   `Bagge2021_visc.nc` has `log10_visc` on `depth`, with 8 structures along
   `structure`.

## Design

### Two levels: run and domain

FastIsostasy is regional and per domain. VILMA2 is one global Earth per run,
shared by every domain (two in `yelmox_bipolar`), like `bsl` today.
`kryos_isostasy` therefore has two types:

```fortran
type isos_run_class          ! driver-owned, one per run (replaces bsl_class in drivers)
    character(len=32)  :: method        ! "fastisostasy" | "vilma2"
    type(bsl_class)    :: bsl           ! fastisostasy: unchanged bsl modes
    type(vilma_run)    :: vl            ! vilma2: solid_earth + global load (stub-able)
    real(wp)           :: bsl_now       ! [m] the bsl every domain uses
end type

type isos_domain_class       ! dom%isos, one per domain
    character(len=32)  :: method
    type(isos_class)   :: fi            ! fastisostasy only
    real(wp), allocatable :: dz_bed(:,:), dz_ss(:,:)   ! on grid_isos: what Yelmo gets
end type
```

The time step is split in three, so several domains can feed one global
solve:

```fortran
call isostasy_load(dom, ir, ts)    ! per domain
call isostasy_solve(ir, ts)        ! once per run
call isostasy_land(dom, ir)        ! per domain -> dom%isos%dz_bed/dz_ss
```

| | fastisostasy | vilma2 |
|---|---|---|
| `isostasy_load` | remap H_ice, dzbdt_corr; `isos_update` (as now) | remap the domain's H_ice onto the global grid; splice into the load |
| `isostasy_solve` | no-op (bsl is updated in `isos_update`, as now) | `solid_earth_update(se, H_global, dt)` |
| `isostasy_land` | `dz_bed = w + we`, `dz_ss = z_ss - bsl` | remap `se%dz_bed`, `se%dz_ss - bsl` back to the domain's `grid_isos` |

`couple_isostasy_to_yelmo` then reads only `dom%isos%dz_bed/dz_ss` and
`ir%bsl_now`, with the same formulas as now. Under `vilma2`, the domain's
`grid_isos` defaults to the VILMA2 global grid, so `isostasy_land` is a
plain copy and the coupler does the single global → Yelmo remap.

Phase 1 of the work below makes this split with FastIsostasy only and must
be bit-identical.

### Parameters

Two groups, split by owner:

- `[vilma]`: the full VILMA2 `&vilma` group (physics and numerics), read by
  `vilma_par_load(se%par, path_par, defaults_file="input/vilma_defaults.nml")`.
  The VILMA2 standalone `&ctl` group is not used.
- `[isos_vilma]` (yelmox, in `yelmox_defaults.nml`): coupling and data.

```
[isos_vilma]
    grid          = "GLOBAL-T64"     ! global Gauss grid (maps/grid_<name>.txt)
    z_bed_ref_file = "fesmdata/Global/{grid_global}/{grid_global}_TOPO-<global product>.nc"
    z_bed_ref_name = "z_bed"
    H_ice_ref_name = "H_ice"
    ice_ext       = "ref"            ! "ref" | "file": ice outside the domains
    ice_ext_file  = "fesmdata/Global/{grid_global}/{grid_global}_IceHistory-ICE6G_C.nc"
    ice_ext_name  = "H_ice"
    dt_ice_ext    = 100.0            ! [yr] how often ice_ext_file is read and interpolated
```

The 3D viscosity file stays a `&vilma` parameter (`visc_3d_file`, new
`visc_3d_structure`). In yelmox's `input/vilma_defaults.nml` its default
is FesmData's `Bagge2021_visc.nc`. Whether 3D viscosity is on stays
VILMA2's own default (`l_visc_3d`).

`input/vilma_defaults.nml` is copied from VILMA2's `input/` and kept in
sync (the header says so, as in `chion_defaults.nml`). Only the 3D
viscosity path differs. `scripts/check_par_nml` validates `[vilma]`
against it.

### Global grid

A Gauss grid named for its truncation (`GLOBAL-T64`, …), both as
`maps/grid_<name>.txt` (cdo `gaussian`) and as a FesmData `gaussian` grid
(`shared/domains.toml`). Its `nlat × nphi` must equal VILMA2's SHT grid
(defaults `nlat = 2 lmax + 2`, `nphi = 4 lmax`; 256 × 130 for `lmax = 64`),
so `grid_is_gauss` is true and VILMA2 does no remap of its own. Init checks
this and stops otherwise. FesmData currently gives `nlat` only and takes
`nphi = 2 nlat` (T170 = 512 × 256), so it needs an `nphi` option, or
`&vilma nphi` is set to `2 nlat`.

FesmData global products (topography, ice history) do not need to be on the
Gauss grid: they are read on a `GLOBAL-<res>DEG` lon-lat grid
(`{grid_global}`) and remapped by the coupler.

The remap method choice (`remap_method_smooth`) compares `dx` in Yelmo
units, which means nothing for a lon-lat grid. Under `vilma2`, the load is
remapped `con` (mass) and the fields come back `bilin`, which is VILMA2's
own convention.

### Global load and reference

- **Load:** `H_global = f·H_domain + (1 - f)·H_ext` on the Gauss grid,
  summed over domains. `f` is the domain footprint fraction (a field of ones,
  remapped `con`; confirm that the coupler leaves cells outside the domain
  at 0).
- **`H_ext`:** `ice_ext = "ref"`: the reference ice `h_ice_eq`, held fixed.
  `"file"`: a global ice history (default ICE-6G_C from FesmData
  IceHistory, Global), read every `dt_ice_ext`, interpolated in time and
  remapped to the Gauss grid. Outside its time range it stops with an
  error.
- **Reference** (`z_bed_eq`, `h_ice_eq`): the global present-day
  topography (`z_bed_ref_file`), with each domain's Yelmo reference
  (`z_bed_ref`, `H_ice_ref`) spliced in the same way. Since only anomalies
  reach Yelmo, the reference only sets VILMA2's ocean function and
  coastlines. `check_isostasy_reference` becomes a footprint-masked
  comparison on the Gauss grid.
- **Ice:** VILMA2 applies flotation itself (`ocean_function`), so it gets
  Yelmo's total `H_ice`, as FastIsostasy does.

### Sea level

- `fastisostasy`: `bsl_class` as now (const / file / fastiso).
- `vilma2`: `ir%bsl_now = se%bsl` (the barystatic rise, gap 2). Any
  `[bsl]` settings are ignored, and a log line says so. `z_sl` gets VILMA2's
  full spatial sea-surface change (geoid + rotation), so it is
  gravitationally self-consistent.

### Time

One cadence per run (`[isostasy] dt`), replacing the per-domain `dt_isos`.
`[comps] with_isostasy` stays per domain. The bed is held between solves, as
now. VILMA2 integrates `dt` internally with its adaptive stepper.

### Restart and output

- Run bundle: `bsl_restart.nc` (fastisostasy) or `vilma_restart.nc`
  (vilma2, `vilma_restart_write`). Per-domain `isos_restart.nc` is written
  for fastisostasy only.
- Per domain, `isos.nc` / `isos_ts.nc` always carry the neutral fields
  (`z_bed`, `z_ss`, `dz_bed`, `dz_ss`, `bsl`); `w`, `we` only for
  fastisostasy.
- Run level, vilma2 only: `vilma.nc` (rsl, dz_bed, dz_ss, C, h_ice on the
  Gauss grid) and `vilma_ts.nc` (bsl, worst SLE mass residual, solver time).

### Spin-up

The cold start calls `solid_earth_init` with the global reference, then
`solid_earth_spinup(se, H_global(t0))` when `[vilma] equil_time_max > 0` or
`pre_spinup_1d`. This is VILMA2's own LGM-memory spin-up, set entirely by
`[vilma]`.

### Build

- `vilma=1`: `VILMAROOT` (configme optional dependency, like rembo1),
  `INC_VILMA`/`LIB_VILMA`, plus SHTns from fesm-utils (`make shtns`).
  Compiles `libs/isos_vilma.f90`, the adapter that owns `type(vilma_run)`.
- default: compiles `libs/isos_vilma_stub.f90` with the same interface;
  `method = "vilma2"` stops at init with "rebuild with vilma=1".
- `libvilma.a` is a file target, and host objects depend on it (the
  stack-wide archive-target convention).

## Changes outside yelmox

VILMA2:

1. Export `dz_bed` (= u, incl. rotational uplift) and `dz_ss` (= rsl + u)
   on the host and Gauss grids, from the converged `u_lm`/`N_lm`.
2. Barystatic sea level as the ocean mean of `C·rsl`. Either `se%bsl`
   changes meaning, or a new `se%gmsl` is added and the current diagnostic
   is renamed.
3. Read FesmData `Bagge2021_visc.nc` (`log10_visc`, `depth`, new
   `visc_3d_structure`). Keep or drop the isostasy_data format.
4. `make vilma-static` → `libvilma/include/libvilma.a` + `.mod`, `openmp`
   variants, as in FastIsostasy.
5. Build and tests against fesm-utils dev; update the README.

FesmData:

1. A global present-day topography product. `Topo` has none yet; the
   domains use BedMachine/Bedmap + IceBoost + GEBCO 2026 merges. Proposal:
   the same merge for `Global` (BedMachine Antarctica v4 and Greenland v6,
   IceBoost v2.0, GEBCO 2026), with `z_bed`, `H_ice`, `z_srf`.
2. Products on the Global lon-lat grids: the topography above and
   `IceHistory-ICE6G_C` (exists). `Bagge2021_visc` is on its own Gaussian
   grid, read directly by VILMA2.
3. Optional: `nphi` for `gaussian` grids, if the yelmox grid files are also
   made by FesmData.

## Phases

- **0 (vilma2, FesmData):** the changes above, with standalone tests
  (benchmarks still pass; `rsl = dz_ss - dz_bed` to round-off; 3D run with
  the FesmData file matches the isostasy_data file for structure 6).
- **1 (yelmox):** add `kryos_isostasy` with FastIsostasy only, move `bsl`
  under `isos_run_class`, and switch the drivers to load/solve/land.
  Bit-identical on the Levante test cases.
- **2 (yelmox):** VILMA2 adapter + stub, `vilma=1` build,
  `input/vilma_defaults.nml`, `[isos_vilma]`, Gauss grid file, global load
  splice (`ice_ext` ref/file), restart, output.
- **3 (tests):** VILMA2 standalone vs yelmox with prescribed ice (same
  load); ANT deglaciation fastisostasy vs vilma2; bipolar (two domains, one
  Earth); exact restart.

## Open questions

1. Default `visc_3d_structure`: 6 (`v_1.0_s16`, identical to the current
   isostasy_data file) or 3 (`v_0.4_s16`, the reference structure in Bagge
   et al. 2021)?
2. Default `lmax` / global grid for yelmox (VILMA2 default 64 →
   `GLOBAL-T64`)?
