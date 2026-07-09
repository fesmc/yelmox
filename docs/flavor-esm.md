---
title: "yelmox_esm"
---

A single-domain driver in which **Earth-System-Model (ESM) output replaces
snapclim/smbpal** as the climate and ocean forcing. The driver owns an
`esm_forcing_class` (from `libs/esm.f90`) instead of a snapclim object, runs it as
a first-class component on its own `grid_clim`, and remaps its fields onto the
consumer grids just like any other module.

- **Program:** `yelmox_esm/yelmox_esm.f90` + `yelmox_esm/yelmox_esm_output.f90` + `libs/yelmox_domain.f90`.
- **Build:** `make yelmox_esm`
- **Config:** `yelmox_esm/yelmox_esm_Antarctica.nml` (and `..._nudge.nml`).

## What's distinct

- **ESM forcing, not snapclim.** `domain_init` is called with
  `init_climate = .false.`, which skips `snapclim_init`; the driver instead calls
  `esm_forcing_init(esm, ..., grid_clim, ...)`. The `esm_forcing_class` holds the
  ESM atmosphere/ocean fields (tas, pr, ocean temperature, and derived inputs).
- **Calendar timeline.** The timeline is initialized with `cal = .true.` and
  `time_ref = 2000.0`, i.e. absolute calendar time (`const_cal`), the ESM
  convention — as opposed to the relative (`time_bp`) timeline used by the
  paleo-oriented flavors.
- **ESM-specific couplers.** Besides the shared `couple_smb_to_yelmo` /
  `couple_marine_to_yelmo`, the driver lands an ESM-owned Yelmo input (a discharge
  field `Qd`) via `couple_esm_extras_to_yelmo`, just before the ice sheet runs.
- **Dedicated climate/marine steps.** `step_climate_esm` and
  `step_marine_shelf_esm` replace the generic `step_climate` / `step_marine_shelf`,
  reading from the `esm_forcing_class` rather than snapclim.

Configuration is read by `esm_ctl_load` from the `[esm]` group: `par_file`,
`experiment`, `esm_name`, and switches `use_esm` / `use_smb` / `use_var` /
`use_proj` / `use_hist`, plus physical parameters (`lapse`, `f_p`, `f_ocn`,
`f_polar`, `dT_threshold`, `grid_src`).

## Stepping order

Main loop (per timestep):

```fortran
call bsl_update(bsl, ts%time_rel)          ! shared sea level, once

call step_optimize(dom, ts)
call step_isostasy(dom, ts, bsl)
call couple_esm_extras_to_yelmo(dom, esm)  ! land ESM-owned Yelmo input (Qd)
call step_icesheet(dom, ts)                ! shared couplers + yelmo_update
call refresh_htopo(dom)

call step_climate_esm(dom, esm, ec, ts)    ! ESM climate + SMB (replaces step_climate)
call step_marine_shelf_esm(dom, esm, ec, ts)
```

This mirrors the single-domain [`yelmox`](flavor-yelmox.md) sequence, with the
generic climate/marine steps swapped for their `*_esm` equivalents and the extra
`couple_esm_extras_to_yelmo` call placed just before `step_icesheet` so the ESM
discharge input is present when Yelmo updates. Cold start builds the initial state
through `esm_cold_start`; restart uses the shared `domain_startup` followed by one
`step_climate_esm` / `step_marine_shelf_esm` to prime the forcing.

## Forcing

Climate transience comes from the ESM data stream itself (historical/projection
experiments selected in `[esm]`), so the `tsgen` `[tsforcing]` mechanism from the
single-domain [`yelmox`](flavor-yelmox.md) is **not** used here.
