---
title: "ISMIP7 projections"
---

Greenland and Antarctic projections for ISMIP7 run in two steps: a present-day
**optimization spin-up**, which fits the basal friction (`cb_ref`) and an ocean
thermal-forcing correction (`tf_corr`) to the observed ice thickness and writes a
restart bundle, then the **scenarios**, which start from that bundle in 2015 and
run under the ESM forcing (historical + projection) to 2300.

- **Program:** `yelmox` (`make yelmox`) with the ESM climate backend
  (`[comps] climate = "esm"`, see [ESM forcing](climate-esm.md)).
- **Config:** `yelmox/yelmox_esm_Greenland.nml`, `yelmox/yelmox_esm_Antarctica_ismip7.nml`.
- **Forcing data:** `input/esm/esm_grl_ismip7.nml`, `input/esm/esm_ant_ismip7.nml`
  (see [Database namelists](database-namelists.md)).
- **Scripts:** `scripts/ismip7/` -- workflows `run_ismip7_greenland.sh`,
  `run_ismip7_antarctica.sh`; tuned spin-ups `opt_grl.sh`, `opt_grl_ismip.sh`,
  `opt_ant.sh`.

The par files on their own run a 10-yr forcing-only smoke test (`[sim] opt =
False`, `[spinup] time_end = 10`); the scripts set the run phase and length.

## Running the workflow

Run from the yelmox root, one step after the other: the scenarios read the
spin-up's last restart bundle.

```bash
scripts/ismip7/run_ismip7_greenland.sh spinup
```

```bash
scripts/ismip7/run_ismip7_greenland.sh scenarios
```

`STAGE=1` creates the run folders and the SLURM script without submitting them
(a dry run). The settings are at the top of each script:

| Setting | Greenland | Antarctica | |
|---|---|---|---|
| `GRID` | `GRL-8KM` | `ANT-16KM` | `[domain] grid_hub` |
| `GCM` | `CESM2-WACCM` | `CESM2-WACCM` | `[esm] esm_name` |
| `SCENARIOS` | `ssp126 ssp370 ssp585` | `ssp585` | `[esm] experiment` |
| `SPINUP_YEARS` | 15000 | 20000 | `[spinup] time_end`, matches `[opt] cf/tf_time_end` |
| `PROJ_END` | 2300 | 2300 | `[transient] time_end` |
| `OUTROOT` | `output/ismip7_grl` | `output/ismip7_ant` | spin-up in `spinup/`, scenarios in `<ssp>/` |

The spin-up passes `ctrl.run_step=spinup sim.relax=True sim.opt=True`; each
scenario passes `ctrl.run_step=transient`, `esm.use_esm/use_hist/use_proj=True`
and `ctrl.restart` = the bundle `<OUTROOT>/spinup/restart-<SPINUP_YEARS/1e3>-kyr`
as an absolute path (the run starts inside its own folder).

**Other grids.** The full input set exists for `GRL-8KM`, `GRL-16KM` and
`ANT-8KM`, `ANT-16KM`, `ANT-32KM` (GRL-4KM: inputs present since 2026-10, not yet
tested). The Antarctic par file reads the Earth structure of ANT-16KM
(`[isos] rheology_file`): with another grid, also pass
`isos.rheology_file=isostasy_data/earth_structure/yelmo/<GRID>_GIA_HR24.nc`
(as `opt_ant.sh` does).

**Short spin-ups.** The spin-up relaxes the ice thickness towards the
observations at first (`[relax] tau1/2`, until `time2`: 6 kyr for
Greenland). A spin-up that ends before then hands the scenarios ice that the
relaxation was holding in place, e.g. thin margin ice under negative SMB, which
is then lost within the first years.

## Tuned spin-ups

`opt_grl.sh`, `opt_grl_ismip.sh` (L. Gutierrez Gonzalez) and `opt_ant.sh` are
15-kyr optimization spin-ups with tuned parameters on top of the par files:
the SMB taken from the ESM (`comps.surface_method=climate`), `cb_ref` started from
the till friction of the bed (`opt.cf_init=-1`, the `ytill` parameters), a short
relaxation (`relax.time1/2 = 100`), equilibrium calving, DIVA dynamics and
shelf enhancement 0.5; `opt_grl_ismip.sh` uses von Mises calving (`vm-l19`) and a
linear melt law. To run the scenarios from one of them, set `BUNDLE` in the
workflow script to its last restart bundle.

## Forcing

The atmosphere and ocean are a reference climatology plus ESM anomalies relative
to `time_esm_ref` (1960--1989), applied over `time_hist` (1850--2014) and
`time_proj` (2015--2300).

| | Greenland | Antarctica |
|---|---|---|
| Reference atmosphere | MAR v3.11 (ERA), 1961--1990 monthly | RACMO2.3p2 (ERA5), 1985--2014 (`time_ref`) |
| Reference ocean | ORAS4 1981--2010 annual mean | meltMIP observational climatology (extrapolated) |
| ESM atmosphere | ISMIP7 SDBN1 `tas`, `pr`, `acabf`, `dacabfdz` | same, `v2` |
| ESM ocean | ISMIP7 `tf`, `so` (2D) | ISMIP7 `thetao`, `so` (3D, `v3`) |
| Subglacial discharge | ISMIP7 `sgd` | -- |
| Shelf melt | `bmb_method = "anom"`, `tf_method = 2` | `bmb_method = "quad-nl"`, `tf_method = 1` |
| Topography, basins | ISMIP7 obs v1.3 (GrIMP surface) | ISMIP7 obs v1.1 (Bedmap3) |

ESMs and scenarios on Levante (`ice_data/ISMIP7/<domain>/<grid>/`):

| ESM | Greenland | Antarctica |
|---|---|---|
| CESM2-WACCM | ctrl, historical, ssp126, ssp370, ssp534-over, ssp585 | historical, ssp585 |
| MRI-ESM2-0 | historical, ssp126, ssp370, ssp585 | historical, ssp126, ssp370, ssp585 |

## Output

The shared per-module files (`yelmo.nc`, `yelmo_ts.nc`, `isos.nc`, `mshlf.nc`,
`smbpal.nc`, `htopo.nc`; Antarctica also `yelmo_ts_APIS/WAIS/EAIS.nc`) and the
restart bundles `restart-<kyr>-kyr/`. For ISMIP7 submissions, `[output]
write_cmip = True` writes `yelmo_cmip.nc` and `yelmo_ts_cmip.nc` every `dt_cmip`;
`write_clim = True` adds the ESM forcing (`esm.nc`, `esm_ts.nc`).
