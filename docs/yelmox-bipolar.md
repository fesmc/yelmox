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
  `active_south`); with neither, the ocean box model runs alone (`active_obm`,
  forced by `atm2obm` or not; `ism2obm` and `ocn_obm<sfx>` off).
- **Shared `bsl`** — one barystatic sea level for the run, restored once at startup
  (`bsl_startup`) and written to the run restart bundle.
- **Shared OBM** (`obm`) — an ocean box model stepped once per timestep, configured
  via `obm_ctl_load` and coupled to the domains through `obm_coupling.f90`
  (`obm_coupling_init`, `obm_exchange`). The OBM writes its own 1D output; its restart
  (`obm_restart.nc`) goes into the run restart bundle next to `bsl_restart.nc`.
- **One restart bundle per run** — `restart-<kyr>-kyr/` holds the shared
  `bsl_restart.nc` and `obm_restart.nc`, and each domain in a subfolder named
  after it (`Greenland/`, `Antarctica/`). `[ctrl] restart` is the bundle to start
  from.

Climate/SMB per domain is **snapshots + chion**, as in the single-domain driver;
the OBM's contribution is folded into the ocean forcing — `obm_exchange` applies
each domain's `[ocn_obm_<sfx>]` (the OBM ocean temperature into `to_ann`) before
the marine-shelf step reads it.

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

call obm_exchange(oc, obox, idx, dom_north, dom_south, ...)  ! atm->obm, ism->obm freshwater,
                                                            ! hysteresis forcing, obm->ism (ocn_obm)

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

One global index drives the box model and, with regional gains, the domains
(`yelmox_bipolar_Bipolar.nml`):

- **Indices** (`[indices]`, `[idx_<name>]`) are shared by both domains and the
  OBM: advanced once per step, written to each domain's 1D file, restarted from
  the run bundle (`idx_<name>_restart.nc`), with the kill switch and
  `restart_every_df` bundles as in `yelmox`. Here `at` (glacial index, the
  north snapshot axis) and `gmt` (global mean air temperature anomaly [K],
  `const` 0; a series instead for a transient).
- **`atm2obm`** (`[ctrl]`): the box air temperatures from `[obm_forcing]`,
  `theta_i = theta_i_init + ampl_i × gain × index` for north, tropics and south
  (`ampl_i`: `[nautilus] thermal_ampl_*`, Zickfeld et al., 2004) and the vapour
  fluxes `phin`, `phit` linear in it. It needs no ice-sheet domain; the index
  must be in the registry (checked at init).
- **Domains**: `[anom_<sfx>]` on the same index with `f_ta` the regional
  amplification (north 1.07, south 0.86): explicit in the par file, zero while
  `gmt` = 0. North: snapshots pd/picontrol/lgm on the axis `at`, ocean from
  `[ocn_fraction_north]`; south: present day only.
- **`obm2ism`**: per domain, the climate modifier `[ocn_obm_<sfx>]`
  (see [`yelmox`](yelmox.md)): `to_ann` = the box temperature (`absolute`) or
  the reference + the box anomaly (`anomaly`), applied by `obm_exchange` after
  both climate steps, and once at the start (`obm_ocean_init`, after the box
  model's cold start or restart, then the shelf step again), so the first
  record and the first step see the box. It needs `active_obm`, and replaces `[ocn_fraction_<sfx>]`
  (refused together).
- **`ism2obm`**: freshwater flux of each domain into the boxes (`fn`, `fs`).

`obm_restart.nc` (run bundle) holds the box state, the fluxes and the box air
temperatures; the `*_init` values stay those of the par file, the reference of
`atm2obm` and of the `anomaly` mode.
