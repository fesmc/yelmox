---
title: "REMBO climate"
aliases:
  - flavor-rembo.html
---

Runs with **REMBOv1**, an energy/moisture-balance regional atmosphere with an
integrated surface-mass-balance scheme, use the `yelmox` program with the REMBO
climate backend: `[comps] climate = "rembo"`. REMBO supplies the atmosphere and
the surface mass balance; the ocean comes from the snapshot backend.

- **Program:** `yelmox/yelmox.f90`, backend in `libs/yelmox_climate.f90`, with the
  REMBO adapter `libs/climate_rembo.f90`.
- **Build:** `make yelmox rembo=1` (links rembo1; prerequisite `rembo-static`).
  Without `rembo=1`, a stub (`libs/climate_rembo_stub.f90`) stops the run if
  `climate = "rembo"`.
- **Config:** `yelmox/yelmox_rembo_Greenland.nml`, plus REMBO's own parameters in
  `yelmox/rembo_Greenland.nml` (staged into the run folder via `.runme/info.json`).
- **Script:** `scripts/rembo/run_rembo.sh [output_dir]` (default `output/rembo`;
  runs locally, `runopts='-rs ...'` submits to the queue).

## Configuration

- **Grid.** REMBO runs on the grid it was compiled for (Greenland, GRL-16KM),
  which must be `grid_clim`; the backend checks this at start-up.
- **Surface mass balance.** `[comps] surface_method = "climate"`: REMBO's smb and
  surface temperature, at the current surface; any other `surface_method` stops
  the run at start-up (REMBO gives annual fields only).
- **Update cadence.** `[comps] dt_clim = dtt` (10 yr): REMBO updates its energy
  balance and its smb on its own intervals (`dtime_emb = 100`, `dtime_smb = 1`
  in `rembo_Greenland.nml`), with the insolation at the forcing time
  `time_rel`.
- **Ocean.** The snapshot backend (see [Snapshots](climate-snapshots.md))
  with an ocean-only reference snapshot: `[forcing] reference = "pd"`, the
  `pd` snapshot binding the ORAS5 1981-2010 climatology (`oras5_1981_2010_to`,
  `oras5_1981_2010_so` of `input/clim_db.nml`). The reference must bind `to`
  and no atmosphere (checked at start-up); `so` is optional (salinity 0
  without it). The ocean takes the ocean terms of `[forcing]` as for the
  snapshots backend (a uniform term on `to`, e.g. `coef = "hyst"`,
  `f_to = 0.25`).
- **Transient forcing.** The uniform terms on `tas` (`[term_<name>]
  pattern = "uniform"`) are an input of REMBO, not applied after it: their
  sum `f_ta·f` is REMBO's summer air temperature anomaly `dT_sum`; the annual
  mean anomaly REMBO applies with it, `dT_ann = (1 + T_wintfac)/2·dT_sum`
  (its seasonal cycle, `T_wintfac` of `rembo_Greenland.nml`), is written as a
  diagnostic. `pr` in their fields is ignored, with a notice at start-up
  (REMBO computes the precipitation); a uniform term on `pr` alone is refused
  (no effect). Snapshot terms on the atmosphere are refused. Off in the par
  file (`[forcing] coefs = ""`, `terms = ""`).
- **Transforms.** `[taper]` and `[bias]` are refused at start-up (they change
  the atmosphere after the backend, and REMBO computes its own smb from
  annual fields); so is a derived `air_anomaly` coefficient (no sea-level air
  temperature).
- **Cold start.** REMBO is equilibrated (10 years) before its first update.

## Output

The shared per-module files (`yelmo.nc`, `isos.nc`, `mshlf.nc`, ...), plus, with
`[output] write_clim`, `rembo.nc` (annual and summer air temperature,
precipitation, smb, surface temperature, ocean profiles) and `rembo_ts.nc` (the applied anomalies, the smb integrated
over the ice sheet and the accumulation-area ratio). The coefficients go to
the Yelmo 1D file (`coef_<name>`), as for the other climates. The library
holds the hysteresis setup, `coef_hyst` (PI42 with kill switch) and the
uniform `term_hyst` on `tas` and `pr`: `-p forcing.coefs=hyst
forcing.terms=hyst` runs it with rembo (its `tas` part is REMBO's `dT_sum`,
its `pr` part ignored; a second uniform term on `to` for the ocean).
Restart bundles hold REMBO's
`rembo_restart.nc`; REMBO reads its restart as set in `rembo_Greenland.nml`.
