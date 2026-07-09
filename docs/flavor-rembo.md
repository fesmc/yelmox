---
title: "yelmox_rembo"
---

A single-domain (Greenland) driver in which **REMBOv1** — an energy/moisture-balance
regional atmosphere with an integrated surface-mass-balance scheme — **replaces the
snapclim + smbpal atmosphere/SMB**. The ocean forcing still comes from snapclim.
REMBO's SMB is staged into the shared SMB carrier so the generic
`couple_smb_to_yelmo` lands it on the Yelmo grid like any other SMB module.

- **Program:** `yelmox_rembo/yelmox_rembo.f90` + `yelmox_rembo/yelmox_rembo_output.f90` + `libs/yelmox_domain.f90`.
- **Build:** `make yelmox_rembo` (links the REMBO stack; prereq `rembo-static`).
- **Configs:** `yelmox_rembo/yelmox_rembo_Greenland.nml` (main) + `yelmox_rembo/rembo_Greenland.nml` (REMBO's own parameters, staged into the run dir).

## What's distinct

- **REMBO atmosphere + SMB.** REMBO is driver-owned (module-global `rembo_ann`),
  initialized with `rembo_init`, and advanced by `step_rembo`, which calls
  `rembo_update` (or `rembo_equilibrate` on the first cold-start call). REMBO runs
  on `grid_clim` (= the Yelmo grid in the single-grid Greenland setup) and works
  internally in double precision.
- **SMB via the shared carrier.** `step_rembo` copies `rembo_ann%smb` / `%T_srf`
  into `dom%smb%ann%smb` / `%tsrf`, so the standard `couple_smb_to_yelmo` (inside
  `step_icesheet`) remaps it conservatively to Yelmo with the water-equiv→ice-equiv
  scaling and the optional `lim_pd_ice` limiter — no REMBO-specific coupler needed.
- **Ocean still from snapclim.** `step_rembo` also calls `snapclim_update` on
  `grid_clim` for the ocean forcing (optionally adding the hysteresis ocean anomaly
  when `ocn_type = "const"`). Atmosphere/SMB no longer come from snapclim/smbpal.

## Stepping order

Main loop (per timestep):

```fortran
call update_hyster_forcing()     ! dT_summer / dT_ann / dT_ocn from the hyster module
call bsl_update(bsl, ts%time_rel)

call step_optimize(dom, ts)
call step_isostasy(dom, ts, bsl)
call step_icesheet(dom, ts)      ! couplers (smb/isos/marine) + yelmo_update
call refresh_htopo(dom)
call step_rembo()                ! REMBO atmosphere/SMB (-> dom%smb%ann) + snapclim ocean
call step_marine_shelf(dom, ts)
```

Like the other multigrid flavors, this inlines the `step_*` primitives and
substitutes `step_rembo` for the generic `step_climate`. REMBO and marine run
*after* `step_icesheet`/`yelmo_update`, so their output is consumed on the next
step (the standard one-step coupling lag). During a hysteresis ramp the driver also
shrinks the main timestep and REMBO's internal `dtime_emb` for stability.

## Forcing

::: {.callout-note}
## REMBO still uses the legacy `hyster` module
Transient forcing here comes from the older `hyster` module (`hyster_init`,
`hyster_calc_forcing` → `dT_summer` / `dT_ann` / `dT_ocn`), **not** the `tsgen`
`[tsforcing]` mechanism used by the single-domain [`yelmox`](flavor-yelmox.md).
Porting REMBO's forcing to `tsgen` is future work; for now it reads the `[hyster]`
namelist group and `ctrl.use_hyster` / `f_ta` / `f_to`.
:::

## Output

REMBO runs write `yelmo2D.nc` (heavy 2D) and `yelmo-rembo.nc` (small 1D + 2D,
including REMBO's `T_ann` / `T_jja` / `pr` diagnostics), plus a restart bundle with
a `rembo_restart.nc`.
