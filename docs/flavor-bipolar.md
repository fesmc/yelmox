---
title: "yelmox_bipolar"
---

A **two-domain** driver that runs a northern and a southern ice-sheet domain
together, coupled through a **shared barystatic sea level** and a **shared Ocean
Box Model (OBM)** that exchanges freshwater flux and ocean temperature between the
hemispheres. Each domain is a full `ice_domain` (the same one the single-domain
[`yelmox`](flavor-yelmox.md) uses); the driver interleaves the `step_*` primitives
across both domains plus the OBM.

- **Program:** `yelmox_bipolar/yelmox_bipolar.f90` + `yelmox_bipolar/obm_coupling.f90` + `libs/yelmox_domain.f90`.
- **Build:** `make yelmox_bipolar` (links the OBM stack, `$(obm_libs)`).
- **Config:** `yelmox_bipolar/yelmox_bipolar_Bipolar.nml`.

## What's distinct

- **Two `ice_domain`s** (`dom_north`, `dom_south`), each set up via `setup_domain`
  → `domain_startup`. Either can be individually deactivated (`active_north` /
  `active_south`).
- **Shared `bsl`** — one barystatic sea level for the run, restored once at startup
  (`bsl_startup`) and written to the run-root restart bundle.
- **Shared OBM** (`obm`) — an ocean box model stepped once per timestep, configured
  via `obm_ctl_load` and coupled to the domains through `obm_coupling.f90`
  (`obm_masks_init`, `obm_exchange`). The OBM writes its own 1D output and restart.

Climate/SMB per domain is still **snapclim + smbpal**, exactly as in the
single-domain driver; the OBM's contribution is folded into the ocean forcing —
`obm_exchange` writes the OBM ocean temperature back into each domain's snapclim
`to_ann` before the marine-shelf step reads it.

## Stepping order

Main loop (per timestep):

```fortran
call bsl_update(bsl, ts%time_rel)              ! shared sea level, once

if (active_north) call advance_isostasy(dom_north)   ! step_optimize + step_isostasy
if (active_south) call advance_isostasy(dom_south)

if (oc%active_obm) call obm_update(obm, dtt, oc%obm_name)   ! ocean box model, one step

if (active_north) call advance_dynamics(dom_north)   ! step_icesheet + refresh_htopo + step_climate
if (active_south) call advance_dynamics(dom_south)

call obm_exchange(oc, obm, dom_north, dom_south, ...)  ! atm->obm, ism->obm freshwater,
                                                       ! hysteresis forcing, obm->ism ocean temp

if (active_north) call step_marine_shelf(dom_north, ts)  ! reads the obm-updated to_ann
if (active_south) call step_marine_shelf(dom_south, ts)
```

Key ordering points:

- **Isostasy for both domains runs before the OBM step**, so the OBM sees a
  consistent geometry.
- **`obm_update` uses the previous step's** atmospheric/freshwater forcing (a
  one-step lag), then `obm_exchange` distributes the fresh OBM state back to the
  domains before the marine-shelf melt is computed.
- `advance_dynamics` = `step_icesheet` → `refresh_htopo` → `step_climate` (the same
  three primitives the single-domain `yelmox_step` runs, minus optimize/isostasy
  which happen in `advance_isostasy` earlier).

## Forcing

Transient forcing is handled through the OBM/hysteresis machinery
(`obm_exchange`), **not** the `tsgen` `[tsforcing]` mechanism — the driver-owned
`tsgen` forcing currently lives only in the single-domain [`yelmox`](flavor-yelmox.md).
