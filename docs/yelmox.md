---
title: "yelmox"
aliases:
  - flavor-yelmox.html
  - flavors.html
---

The YelmoX program: one ice-sheet domain with FastIsostasy bedrock, a shared
barystatic sea level, a climate backend chosen at runtime, a surface mass balance
and the marine shelf, built on the multigrid `kryos_domain` (see
[Multigrid coupling](multigrid.md)). Its variant
[`yelmox_bipolar`](yelmox-bipolar.md) runs two domains, north and south, coupled
through the sea level and an ocean box model.

| Program | Build | Domains | Climate | Ocean |
|---|---|---|---|---|
| `yelmox` | `make yelmox` (`rembo=1` for REMBO) | one | any backend (`[comps] climate`) | the climate (depth profiles or shelf base) |
| [`yelmox_bipolar`](yelmox-bipolar.md) | `make yelmox_bipolar` | north + south | snapclim (×2) | snapclim + shared ocean box model |

- **Program:** `yelmox/yelmox.f90` (thin driver) + `libs/kryos*.f90` (domain, coupling, startup, output).
- **Build:** `make yelmox` (`make yelmox rembo=1` to link REMBO).
- **Configs:** `yelmox/yelmox_<domain>.nml` (Antarctica, Greenland, North, LIS,
  Pyrenees, SRG, plus `pd_` present-day and paleo variants), `yelmox_esm_*.nml`
  ([ESM forcing](climate-esm.md)), `yelmox_rembo_Greenland.nml`
  ([REMBO](climate-rembo.md)), `yelmox_Greenland_snapesm*.nml`
  ([snapesm](climate-snap.md)).

## Components

| Role | Module | Grid |
|---|---|---|
| Ice sheet | Yelmo | `grid_ice` |
| Isostasy + sea level | FastIsostasy (`isos`) + shared `bsl` | `grid_isos` |
| Climate (atmosphere + ocean) | backend of `[comps] climate`: [snapclim, snapesm](climate-snap.md), [esm](climate-esm.md) or [rembo](climate-rembo.md) | `grid_clim` |
| Surface mass balance | chion (default), smbpal, `smb_simple` or the climate's own | `grid_surface` |
| Sub-shelf melt | marine_shelf | `grid_shelf` |
| Geometry hub | htopo | `grid_hub` (hi-res) |

Each module runs on its own grid, set in `[domain]`; the coupler remaps fields between
grids at the moment of coupling. See [Multigrid coupling](multigrid.md).

### Climate backends

The forcing of a run is set at runtime by its climate backend, `[comps] climate`:

| Backend | Supplies | Page |
|---|---|---|
| `snapclim` | atmosphere + ocean profiles, from climate snapshots blended by indices | [snapclim and snapesm](climate-snap.md) |
| `snapesm` | as snapclim, configured as one blend model over any number of snapshots | [snapclim and snapesm](climate-snap.md) |
| `esm` | reference climatology + Earth-system-model anomalies; atmosphere or smb, ocean at the shelf base, subglacial discharge | [ESM forcing](climate-esm.md) |
| `rembo` | REMBOv1 atmosphere + smb; ocean from snapclim | [REMBO climate](climate-rembo.md) |

## Stepping order

The driver owns the timeline (`ts`) and the shared sea level (`bsl`), and advances
the domain once per step with the coupling sequence written out in the time loop:

```fortran
call step_relax(dom, ts)          ! topography relaxation (relax)
call step_optimize(dom, ts)       ! cb_ref/tf_corr optimization (opt)
call step_isostasy(dom, ts, bsl)  ! bedrock + sea level, this step
call couple_to_yelmo(dom)         ! bedrock now; smb + shelf melt lag one step
call step_icesheet(dom, ts)       ! yelmo_update
call couple_yelmo_to_htopo(dom)   ! hi-res geometry from the models
call step_climate(dom, ts, tsf)   ! climate (dt_clim cadence)
call step_surface(dom, ts)        ! surface mass balance + temperature
call step_shelf(dom, ts)          ! shelf-base melt + temperature
```

`couple_to_yelmo` assembles the Yelmo boundary state: isostasy from this step,
smb and shelf melt from the **previous** step (a one-step coupling lag);
`step_climate`, `step_surface` and `step_shelf` then produce the forcing
consumed on the next step.
The climate is refreshed on the `comps.dt_clim` cadence; the smb every step.

## Transient time-series forcing (`tsgen`)

`yelmox` can drive a spatially-homogeneous, time-varying anomaly into the climate
(atmosphere and/or ocean) from the `tsgen` time-series generator (the modern
replacement for the legacy `hyster` module). It is **driver-owned**: the program
holds a `tsforcing_class` (`tsf`), advances it each step, and passes it to
`step_climate` (and the cold start), which apply its `dTa` / `dTo` / `dSo`.

Two namelist groups control it:

```fortran
&tsforcing
    active = True     ! turn transient forcing on
    f_ta   = 1.0      ! dTa = f_now * f_ta   (atmospheric temperature [K])
    f_to   = 0.0      ! dTo = f_now * f_to   (ocean temperature [K])
    f_so   = 0.0      ! dSo = f_now * f_so   (ocean salinity [psu])
/

&tsgen
    method    = "ramp-time"   ! const | ramp-slope | ramp-time | ramp-time-step | sin | exp | PI42 | ...
    f_min     = 0.0
    f_max     = 5.0
    dt_ramp   = 200.0
    ! ... (see tsgen.f90 for the full parameter set)
/
```

`tsgen` produces a single scalar `f_now`, which the `[tsforcing]` gains map onto
the three anomalies. Time-driven methods (`ramp-*`, `sin`, `const`) are
analytic; feedback methods (`exp`, PI/PID controllers) modulate the forcing rate
from the model response — the response variable passed to `tsgen` is total ice
volume (Gt).

How the anomalies reach the climate depends on the backend:

| Backend | `dTa` / `dTo` / `dSo` |
|---|---|
| snapclim | used only in the `"anom"` methods (`snap.atm_type` for `dTa`, `snap.ocn_type` for `dTo`/`dSo`); ignored in the index-based methods (`snap_1ind_new`, `snap_2ind`, `hybrid`, …), as in the legacy `hyster` contract |
| snapesm | added on top in every configuration |
| esm | not used: the backend has its own forcing |
| rembo | mapped onto REMBO's summer, annual and ocean anomalies (see [REMBO](climate-rembo.md)) |

With `active = False` (the default in the shipped configs) `step_climate`
applies no anomalies.

## Configuration

A run is set by these namelist groups. Every key is required (a missing key stops
the run), except the keys of an unselected `init_method` and the `[relax]` and
`[opt]` groups when their switch in `[sim]` is off.

| Group | Sets |
|---|---|
| `[ctrl]` | the run: its timeline group, calendar years, the restart bundle to start from |
| timeline group (`[ctrl] run_step`) | the timeline: `tstep_method`, `tstep_const`, `time_init`, `time_end`, `dtt` |
| `[domain]` | the domain: name, grid of every component, hub topography and masks (see [Multigrid coupling](multigrid.md)) |
| `[comps]` | the components: which are active, with which model, how often |
| `[sim]` | the conditions of this simulation: cold-start ice state, relaxation, optimization, regional modifications |
| `[relax]`, `[opt]` | the parameters of the relaxation and the optimization switched on in `[sim]` |
| `[output]` | the output files; their intervals are `[tm_1D]`, `[tm_2D]`, `[tm_2Dsm]`, and `[tm_rst]` for restarts |

The models have their own groups (`[yelmo]` and the Yelmo physics groups,
`[isos]`, `[barysealevel]`, `[marine_shelf]`, `[surface_chion]` + `[chion]`, `[smbpal]`, the climate backend's,
...). In `yelmox_bipolar` the domain groups carry the hemisphere suffix
(`[domain_north]`, `[comps_north]`, `[sim_north]`, `[relax_north]`,
`[opt_north]`, `[output_north]`, ...); `[ctrl]` and the timeline are shared.

A minimal skeleton:

```fortran
&ctrl
    run_step        = "ctrl"            ! timeline group ("ctrl" = this group)
    calendar        = False
    calendar_ref    = 1950.0
    restart         = "None"            ! restart bundle folder, or "None" for a cold start
    tstep_method    = "const"
    tstep_const     = 0.0
    time_init       = 0.0
    time_end        = 15e3
    dtt             = 10.0
/

&comps
    with_ice_sheet  = True
    with_isostasy   = True
    with_climate    = True
    with_surface    = True
    with_shelf      = True
    climate         = "snapclim"
    surface_method  = "chion"
    dt_clim         = 10.0
/

&sim
    init_method       = "none"
    init_marine_H     = False
    init_kill_shelves = False
    init_time_thrm    = 0.0
    relax             = True            ! uses [relax]
    opt               = True            ! uses [opt]
    use_negis         = False
    scale_glacial_smb = False
    lim_pd_ice        = False
/
```

### `[ctrl]`

| Key | Values | |
|---|---|---|
| `run_step` | group name | the group holding the timeline: `"ctrl"` for `[ctrl]` itself, or a phase group such as `[spinup]` or `[transient]`, so one par file can hold several run phases |
| `calendar`, `calendar_ref` | bool, [yr CE] | timeline in calendar years, against the reference year `calendar_ref` (ESM runs) |
| `restart` | folder, `"None"` | the restart bundle to start from; `"None"` = cold start |

The timeline group holds `tstep_method` (`"const"`: the time advances and the
forcing time stays at `tstep_const`; `"rel"`: the forcing follows the time, e.g.
for paleo records; `"cal"`: calendar time), `tstep_const`, `time_init`,
`time_end` and the main time step `dtt` [yr].

A restart bundle is a folder `restart-<kyr>-kyr/`, written on the `[tm_rst]`
schedule and always at `time_end`. It holds one restart file per stateful model,
the shared sea level and the transient-forcing state (in `yelmox_bipolar` also the
ocean box model, with each domain in a subfolder named after it). Pass `restart`
as an absolute path: the run starts inside its own folder.

### `[comps]`

The components: which are active, with which model, how often.

| Key | Values | |
|---|---|---|
| `with_ice_sheet`, `with_isostasy`, `with_climate`, `with_surface`, `with_shelf` | bool | components in the coupling sequence: the ice sheet, isostasy, the climate (atmosphere + ocean), the surface (mass balance + temperature) and the shelf base (melt + temperature); `with_surface` and `with_shelf` need `with_climate` |
| `climate` | `snapclim`, `snapesm`, `esm`, `rembo` | the climate backend |
| `surface_method` | `chion` (default), `smbpal`, `smb_simple`, `climate` | chion or smbpal (ITM or PDD), both from the climate's temperature and precipitation; chion runs the ITM or BESSI snowpack (`[chion] model`) in daily steps (`[surface_chion]`, `[chion]`; ITM's parameters are the group `[chion] nml_itm` names: smbpal's `[itm]`, or for Greenland `[itm_chion]` (`alb_ice = 0.31`, calibrated to MAR); BESSI gets the surface shortwave as TOA insolation times `trans_sw`, a constant `wind_speed` and `rel_hum`, and the air pressure from the surface elevation), smbpal also has PDD (the PDD configurations use it); smb_simple (needs a sea-level air temperature: snapclim, snapesm); the climate's own smb (esm, rembo; required by rembo) |
| `dt_clim` | [yr] | climate update interval; `<= 0`: updated only at the cold start |

### `[sim]`

The conditions of this simulation: the cold-start ice state (`init_*`), the
relaxation, the optimization and regional modifications. The `init_*` keys act on
a cold start only; `relax` and `opt` act on a cold start or a restart, with their
times counted from the start of the run (`time_init`).

| Key | Values | |
|---|---|---|
| `init_method` | `none`, `equil`, `recon`, `recon_ref` | cold-start ice state: as initialized; a short equilibration with constant boundaries; the reconstruction `recon_path` as initial ice where `recon_regions` selects; the reconstruction as reference ice only |
| `init_equil_time` | [yr] | `equil`: equilibration time |
| `recon_path`, `recon_var` | path, name | `recon`, `recon_ref`: the reconstruction file (`{domain}`, `{grid_name}` = `grid_ice`) and its ice-thickness variable |
| `recon_regions` | selection expression | `recon`: where its ice is imposed (e.g. `"region:North_America"`) |
| `init_marine_H` | bool | LGM-like marine ice, before `init_method` |
| `init_kill_shelves` | bool | no ice where the present-day bed is ocean |
| `init_time_thrm` | [yr] | then equilibrate with the topography fixed (`0` = off) |
| `relax` | bool | relax the topography towards the reference (`[relax]`) |
| `opt` | bool | optimize the basal friction and the thermal-forcing correction (`[opt]`) |
| `use_negis` | bool | NEGIS basal-friction modification (`[negis]`) |
| `scale_glacial_smb` | bool | reduce negative glacial smb (`[glacial_smb]`) |
| `lim_pd_ice` | bool | extra melt (4 m/yr) outside the present-day ice extent |

A spin-up usually sets `relax` and `opt` together: the relaxation holds the
floating ice and the grounding zone near the observations while the optimization
adjusts the basal friction, and then releases them gradually. A run restarted from
a spin-up sets both to `False` to keep the optimized fields fixed.

### `[relax]`

The relaxation of the topography (`step_relax`). While active, `ytopo.topo_rel`
is `topo_rel` and its timescale ramps from `tau1` to `tau2`; after `time2`, the
`[ytopo]` values of `topo_rel` and `topo_rel_tau` apply again.

| Key | Values | |
|---|---|---|
| `topo_rel` | `ytopo.topo_rel` mode | where the ice relaxes while active (`3`: all points, `4`: the grounding line and grounding zone) |
| `tau1`, `tau2` | [yr] | relaxation timescale until `time1`, and at `time2` |
| `time1`, `time2` | [yr] | end of the `tau1` period, and end of the relaxation |
| `m` | [-] | exponent of the ramp from `tau1` to `tau2` |

### `[opt]`

The optimization (`step_optimize`), towards the observed ice thickness, of the
basal-friction coefficient `cb_ref` (`opt_cf`, between `cf_time_init` and
`cf_time_end`) and of the thermal-forcing correction `tf_corr` of the marine shelf
(`opt_tf`, between `tf_time_init` and `tf_time_end`):

| Key | Values | |
|---|---|---|
| `opt_cf` | `none`, `L21` | `cb_ref` following Lipscomb et al. (2021) |
| `opt_tf` | `none`, `L21`, `L21-points` | one `tf_corr` per basin (`tf_basins`), or one per point (`tf_sigma`, `basin_fill`) |
| `cf_init` | value | initial `cb_ref` on a cold start (`<= 0`: the till friction of the bed, `cb_tgt`) |

The method and its other parameters are described in the Yelmo docs,
[Basal friction optimization](https://fesmc.github.io/yelmo/optimization.html).

### `[output]`

Each module writes its own files, on its own grid, at the `[tm_2D]` (2D),
`[tm_2Dsm]` (small 2D) and `[tm_1D]` (time series) intervals:

| Key | Files |
|---|---|
| `write_yelmo` | `yelmo.nc`, `yelmo_sm.nc`, `yelmo_ts.nc` (and `yelmo_ts_<region>.nc` for the named regions of `[domain]`) |
| `write_isos` | `isos.nc`, `isos_ts.nc` |
| `write_shelf` | `mshlf.nc` |
| `write_surface` | `chion.nc` (chion), `smbpal.nc` (smbpal, climate) |
| `write_clim` | the backend's file: `snap.nc` (snapclim, snapesm), `esm.nc` + `esm_ts.nc`, `rembo.nc` + `rembo_ts.nc` |
| `write_htopo` | `htopo.nc` (the hub) |
| `write_cmip`, `dt_cmip` | `yelmo_cmip.nc`, `yelmo_ts_cmip.nc`, every `dt_cmip` years (marine-shelf fields need `grid_shelf = grid_ice`) |

## Also built from this driver

`make yelmox_glaciers` previously produced a second binary from this same
`yelmox.f90` for mountain-glacier runs; it has been removed. Use `make yelmox`
with the `yelmox_Pyrenees.nml` / `yelmox_SRG.nml` configs instead.
