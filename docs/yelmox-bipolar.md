---
title: "yelmox_bipolar"
aliases:
  - flavor-bipolar.html
---

The two-domain variant of [`yelmox`](yelmox.md): it runs a northern and a southern
ice-sheet domain together, coupled through a **shared barystatic sea level** and a **shared Ocean
Box Model (OBM)** that exchanges freshwater flux and ocean temperature between the
hemispheres. Each domain is a full `kryos_domain`, configured
as in `yelmox` with the hemisphere suffix on its groups (see
[Configuration](yelmox.md#configuration)); the driver interleaves the `step_*` primitives
across both domains plus the OBM.

- **Program:** `yelmox_bipolar/yelmox_bipolar.f90` + `yelmox_bipolar/obm_coupling.f90` + `libs/kryos*.f90`.
- **Build:** `make yelmox_bipolar` (links the OBM stack, `$(obm_libs)`).
- **Config:** `yelmox_bipolar/yelmox_bipolar_Bipolar.nml`.

## What's distinct

- **Two `kryos_domain`s** (`dom_north`, `dom_south`), each set up via `setup_domain`
  → `domain_startup`. Either can be individually deactivated (`active_north` /
  `active_south`).
- **Shared `bsl`** — one barystatic sea level for the run, restored once at startup
  (`bsl_startup`) and written to the run restart bundle.
- **Shared OBM** (`obm`) — an ocean box model stepped once per timestep, configured
  via `obm_ctl_load` and coupled to the domains through `obm_coupling.f90`
  (`obm_masks_init`, `obm_exchange`). The OBM writes its own 1D output; its restart
  (`obm_restart.nc`) goes into the run restart bundle next to `bsl_restart.nc`.
- **One restart bundle per run** — `restart-<kyr>-kyr/` holds the shared
  `bsl_restart.nc` and `obm_restart.nc`, and each domain in a subfolder named
  after it (`Greenland/`, `Antarctica/`). `[ctrl] restart` is the bundle to start
  from.

Climate/SMB per domain is still **snapclim + chion (north) / smbpal PDD (south)**, exactly as in the
single-domain driver; the OBM's contribution is folded into the ocean forcing —
`obm_exchange` writes the OBM ocean temperature back into each domain's snapclim
`to_ann` before the marine-shelf step reads it.

## Stepping order

Main loop (per timestep):

```fortran
call bsl_update(bsl, ts%time_rel)              ! shared sea level, once

if (active_north) then                         ! relaxation + optimization + isostasy, per domain
    call step_relax(dom_north, ts)
    call step_optimize(dom_north, ts)
    call step_isostasy(dom_north, ts, bsl)
end if
(same for dom_south)

if (oc%active_obm) call obm_update(obox, dtt, oc%obm_name)   ! ocean box model, one step

if (active_north) then                         ! ice sheet, hub, climate + smb, per domain
    call couple_to_yelmo(dom_north)
    call step_icesheet(dom_north, ts)
    call couple_yelmo_to_htopo(dom_north)
    call step_climate(dom_north, ts)
    call step_surface(dom_north, ts)
end if
(same for dom_south)

call obm_exchange(oc, obox, dom_north, dom_south, ...)  ! atm->obm, ism->obm freshwater,
                                                       ! hysteresis forcing, obm->ism ocean temp

if (active_north) call step_shelf(dom_north, ts)  ! reads the obm-updated to_ann
if (active_south) call step_shelf(dom_south, ts)
```

Key ordering points:

- **Isostasy for both domains runs before the OBM step**, so the OBM sees a
  consistent geometry.
- **`obm_update` uses the previous step's** atmospheric/freshwater forcing (a
  one-step lag), then `obm_exchange` distributes the fresh OBM state back to the
  domains before the marine-shelf melt is computed.
- The per-domain primitives are the same as in the single-domain `yelmox`; only
  the relaxation + optimization + isostasy part is split off so it runs before the OBM step.

## Forcing

Transient forcing is handled through the OBM/hysteresis machinery
(`obm_exchange`), **not** the `tsgen` `[tsforcing]` mechanism — the driver-owned
`tsgen` forcing currently lives only in the single-domain [`yelmox`](yelmox.md).
