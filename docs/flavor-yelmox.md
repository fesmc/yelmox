---
title: "yelmox (single-domain)"
---

The single-domain program: one ice-sheet domain with FastIsostasy bedrock, a
shared barystatic sea level, a climate backend chosen at runtime, a surface mass
balance and the marine shelf. It is the reference implementation of the multigrid
`kryos_domain`; [`yelmox_bipolar`](flavor-bipolar.md) runs two of them.

- **Program:** `yelmox/yelmox.f90` (thin driver) + `libs/kryos*.f90` (domain, coupling, startup, output).
- **Build:** `make yelmox` (`make yelmox rembo=1` to link REMBO).
- **Configs:** `yelmox/yelmox_<domain>.nml` (Antarctica, Greenland, North, LIS,
  Pyrenees, SRG, plus `pd_` present-day and paleo variants), `yelmox_esm_*.nml`
  ([ESM forcing](flavor-esm.md)), `yelmox_rembo_Greenland.nml`
  ([REMBO](flavor-rembo.md)), `yelmox_Greenland_snapesm*.nml`
  ([snapesm](climate-snap.md)).

## Components

| Role | Module | Grid |
|---|---|---|
| Ice sheet | Yelmo | `grid_ice` |
| Isostasy + sea level | FastIsostasy (`isos`) + shared `bsl` | `grid_isos` |
| Climate (atmosphere + ocean) | backend of `[coupling] climate`: [snapclim, snapesm](climate-snap.md), [esm](flavor-esm.md) or [rembo](flavor-rembo.md) | `grid_clim` |
| Surface mass balance | smbpal, `smb_simple` or the climate's own | `grid_smb` |
| Sub-shelf melt | marine_shelf | `grid_mshlf` |
| Geometry hub | htopo | `grid_hub` (hi-res) |

Each module runs on its own grid, set in `[domain]`; the coupler remaps fields between
grids at the moment of coupling. See [Multigrid coupling](multigrid.md).

## Stepping order

The driver owns the timeline (`ts`) and the shared sea level (`bsl`), and advances
the domain once per step with the coupling sequence written out in the time loop:

```fortran
call step_spinup_tuning(dom, ts)  ! relaxation ramp + cb_ref/tf_corr tuning (opt)
call step_isostasy(dom, ts, bsl)  ! bedrock + sea level, this step
call couple_to_yelmo(dom)         ! bedrock now; smb + shelf melt lag one step
call step_icesheet(dom, ts)       ! yelmo_update
call couple_yelmo_to_htopo(dom)   ! hi-res geometry from the models
call step_climate(dom, ts, tsf)   ! climate (dt_clim cadence)
call step_smb(dom, ts)            ! surface mass balance
call step_marine_shelf(dom, ts)   ! shelf melt
```

`couple_to_yelmo` assembles the Yelmo boundary state: isostasy from this step,
smb and shelf melt from the **previous** step (a one-step coupling lag);
`step_climate`, `step_smb` and `step_marine_shelf` then produce the forcing
consumed on the next step.
The climate is refreshed on the `coupling.dt_clim` cadence; the smb every step.

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
| rembo | mapped onto REMBO's summer, annual and ocean anomalies (see [REMBO](flavor-rembo.md)) |

With `active = False` (the default in the shipped configs) `step_climate`
applies no anomalies.

## Configuration

Besides `[domain]` (the domain definition, see [Multigrid coupling](multigrid.md)),
a run is set by `[coupling]` and `[output]`. In `yelmox_bipolar` each group
carries the hemisphere suffix (`[coupling_north]`, ...). Every key is required (a
missing key stops the run), except the keys of an unselected `init_method`.

### `[coupling]`

| Key | Values | |
|---|---|---|
| `with_ice_sheet`, `with_isostasy`, `with_climate`, `with_marine_shelf` | bool | components in the coupling sequence; `with_climate` covers climate and smb |
| `climate` | `snapclim`, `snapesm`, `esm`, `rembo` | the climate backend |
| `smb_method` | `smbpal`, `smb_simple`, `climate` | smbpal (from the climate's temperature and precipitation); smb_simple (needs a sea-level air temperature: snapclim, snapesm); the climate's own smb (esm, rembo; required by rembo) |
| `dt_clim` | [yr] | climate update interval; `<= 0`: updated only at the cold start |
| `equil_method` | `none`, `opt` | `opt`: spin-up optimization of the basal friction and the thermal-forcing correction (`[opt]`) |
| `restart` | folder, `"None"` | restart bundle to start from; `"None"` = cold start |
| `init_method` | `none`, `equil`, `recon`, `recon_ref` | cold-start ice state: as initialized; a short equilibration with constant boundaries; the reconstruction `recon_path` as initial ice on the `recon_codes` regions; the reconstruction as reference ice only |
| `init_equil_time` | [yr] | `equil`: equilibration time |
| `recon_path`, `recon_var` | path, name | `recon`, `recon_ref`: the reconstruction file (`{domain}`, `{grid_name}` = `grid_ice`) and its ice-thickness variable |
| `recon_codes` | codes | `recon`: the regions where its ice is imposed |
| `init_marine_H` | bool | cold start: LGM-like marine ice, before `init_method` |
| `kill_shelves` | bool | cold start: no ice where the present-day bed is ocean |
| `time_equil_thrm` | [yr] | cold start: then equilibrate with topography fixed (`0` = off) |
| `scale_glacial_smb` | bool | reduce negative glacial smb (`[glacial_smb]`) |
| `lim_pd_ice` | bool | extra melt (4 m/yr) outside the present-day ice extent |
| `use_negis` | bool | NEGIS basal-friction modification (`[negis]`) |

### `[output]`

Each module writes its own files, on its own grid, at the `[tm_2D]` (2D),
`[tm_2Dsm]` (small 2D) and `[tm_1D]` (time series) intervals:

| Key | Files |
|---|---|
| `write_yelmo` | `yelmo.nc`, `yelmo_sm.nc`, `yelmo_ts.nc` (and `yelmo_ts_<region>.nc` for the named regions of `[domain]`) |
| `write_isos` | `isos.nc`, `isos_ts.nc` |
| `write_mshlf` | `mshlf.nc` |
| `write_smb` | `smbpal.nc` |
| `write_clim` | the backend's file: `snap.nc` (snapclim, snapesm), `esm.nc` + `esm_ts.nc`, `rembo.nc` + `rembo_ts.nc` |
| `write_htopo` | `htopo.nc` (the hub) |
| `write_cmip`, `dt_cmip` | `yelmo_cmip.nc`, `yelmo_ts_cmip.nc`, every `dt_cmip` years (marine-shelf fields need `grid_mshlf = grid_ice`) |

Restart bundles follow `[tm_rst]`, plus one at `time_end`.

## Also built from this driver

`make yelmox_glaciers` previously produced a second binary from this same
`yelmox.f90` for mountain-glacier runs; it has been removed. Use `make yelmox`
with the `yelmox_Pyrenees.nml` / `yelmox_SRG.nml` configs instead.
