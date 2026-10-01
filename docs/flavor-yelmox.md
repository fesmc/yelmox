---
title: "yelmox (single-domain)"
---

The canonical driver: a single ice-sheet domain forced by **snapclim** (climate +
ocean) and **smbpal** (surface mass balance), with FastIsostasy bedrock and a
shared barystatic sea level. It is the reference implementation of the multigrid
`kryos_domain` and the template the other flavors specialize.

- **Program:** `yelmox/yelmox.f90` (thin driver) + `libs/kryos*.f90` (domain, coupling, startup, output).
- **Build:** `make yelmox`
- **Configs:** `yelmox/yelmox_<domain>.nml` (Antarctica, Greenland, North, LIS, Pyrenees, SRG, plus `pd_` present-day variants).

## Components

| Role | Module | Grid |
|---|---|---|
| Ice sheet | Yelmo | `grid_yelmo` |
| Isostasy + sea level | FastIsostasy (`isos`) + shared `bsl` | `grid_isos` |
| Climate (atmosphere + ocean) | snapclim | `grid_clim` |
| Surface mass balance | smbpal (or `smb_simple`) | `grid_smb` |
| Sub-shelf melt | marine_shelf | `grid_mshlf` |
| Geometry hub | htopo | `grid_name` (hi-res) |

Each module runs on its own configurable grid; the coupler remaps fields between
grids at the moment of coupling. See [Multigrid coupling](multigrid.md).

## Stepping order

The driver owns the timeline (`ts`) and the shared sea level (`bsl`), and advances
the domain once per step with the coupling sequence written out in the time loop:

```fortran
call step_spinup_tuning(dom, ts)  ! relaxation ramp + cb_ref/tf_corr tuning (opt)
call step_isostasy(dom, ts, bsl)  ! bedrock + sea level, this step
call couple_to_yelmo(dom)         ! bedrock now; smb + shelf melt lag one step
call step_icesheet(dom, ts)       ! yelmo_update
call refresh_hub(dom)             ! hi-res geometry from the models
call step_climate(dom, ts, tsf)   ! climate (dt_clim cadence)
call step_smb(dom, ts)            ! surface mass balance
call step_marine_shelf(dom, ts)   ! shelf melt
```

`couple_to_yelmo` assembles the Yelmo boundary state: isostasy from this step,
smb and shelf melt from the **previous** step (a one-step coupling lag);
`step_climate`, `step_smb` and `step_marine_shelf` then produce the forcing
consumed on the next step.
snapclim is refreshed on the `coupling.dt_clim` cadence; smbpal every step.

## Transient time-series forcing (`tsgen`)

`yelmox` can drive a spatially-homogeneous, time-varying anomaly into snapclim
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
the three snapclim anomalies. Time-driven methods (`ramp-*`, `sin`, `const`) are
analytic; feedback methods (`exp`, PI/PID controllers) modulate the forcing rate
from the model response — the response variable passed to `tsgen` is total ice
volume (Gt).

::: {.callout-important}
## `dTa`/`dTo`/`dSo` require snapclim's `"anom"` mode
The anomalies are only consumed when `snap.atm_type = "anom"` (for `dTa`) and
`snap.ocn_type = "anom"` (for `dTo`/`dSo`). In snapclim's index-based modes
(`snap_1ind_new`, `snap_2ind`, `hybrid`, …) the passed anomalies are **ignored**
— the forcing there comes from the index time-series files instead. This matches
the legacy `hyster` contract.
:::

With `active = False` (the default in the shipped configs) `step_climate`
applies no anomalies and snapclim behaves exactly as before.

## Also built from this driver

`make yelmox_glaciers` previously produced a second binary from this same
`yelmox.f90` for mountain-glacier runs; it has been removed. Use `make yelmox`
with the `yelmox_Pyrenees.nml` / `yelmox_SRG.nml` configs instead.
