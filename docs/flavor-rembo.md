---
title: "REMBO climate"
---

Runs with **REMBOv1**, an energy/moisture-balance regional atmosphere with an
integrated surface-mass-balance scheme, use the `yelmox` program with the REMBO
climate backend: `[comps] climate = "rembo"`. REMBO supplies the atmosphere and
the surface mass balance; the ocean comes from snapclim.

- **Program:** `yelmox/yelmox.f90`, backend in `libs/yelmox_climate.f90`, with the
  REMBO adapter `libs/climate_rembo.f90`.
- **Build:** `make yelmox rembo=1` (links rembo1; prerequisite `rembo-static`).
  Without `rembo=1`, a stub (`libs/climate_rembo_stub.f90`) stops the run if
  `climate = "rembo"`.
- **Config:** `yelmox/yelmox_rembo_Greenland.nml`, plus REMBO's own parameters in
  `yelmox/rembo_Greenland.nml` (staged into the run folder via `.runme/info.json`).
- **Script:** `scripts/rembo/run_rembo.sh`.

## Configuration

- **Grid.** REMBO runs on the grid it was compiled for (Greenland, GRL-16KM),
  which must be `grid_clim`; the backend checks this at start-up.
- **Surface mass balance.** `[comps] surface_method = "climate"`: REMBO's smb and
  surface temperature, at the current surface. REMBO gives annual fields only, so
  smbpal and smb_simple are not available with it.
- **Update cadence.** `[comps] dt_clim = dtt`: REMBO updates its energy balance
  and its smb on its own intervals (`dtime_emb`, `dtime_smb`).
- **Ocean.** The `[snap]` group, as for snapclim.
- **Transient forcing.** `[tsforcing]` maps the tsgen value `f_now` onto REMBO's
  anomalies: the summer air temperature `dT_sum = f_now·f_ta`, the annual
  `dT_ann = 1.3·dT_sum` (REMBO's winter factor 1.6), and the ocean
  `dT_ocn = dT_ann·f_to`, added to a snapclim ocean held at its reference
  (`ocn_type = "const"`).
- **Cold start.** REMBO is equilibrated (10 years) before its first update.

## Output

The shared per-module files (`yelmo.nc`, `isos.nc`, `mshlf.nc`, ...), plus, with
`[output] write_clim`, `rembo.nc` (annual and summer air temperature,
precipitation, smb) and `rembo_ts.nc` (the applied anomalies, the smb integrated
over the ice sheet and the accumulation-area ratio). The tsgen forcing goes to
`yelmo_ts.nc`, as for the other climates. Restart bundles hold REMBO's
`rembo_restart.nc`; REMBO reads its restart as set in `rembo_Greenland.nml`.
