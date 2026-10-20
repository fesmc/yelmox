---
title: "ESM forcing"
aliases:
  - flavor-esm.html
---

Runs forced by **Earth-System-Model (ESM) output** (ISMIP7, TIPMIP, 1pctCO2) use
the `yelmox` program with the ESM climate backend: `[comps] climate = "esm"`.
The backend wraps `libs/esm_forcing.f90`: a reference climatology at the current
surface (lapse rate, precipitation scaling), plus ESM anomalies over historical and
projection periods, optional climate variability and subglacial discharge.

- **Program:** `yelmox/yelmox.f90` (`make yelmox`), backend in `libs/yelmox_climate.f90`.
- **Config:** `yelmox/yelmox_esm_Greenland.nml`, `yelmox/yelmox_esm_Antarctica.nml`.
  The experiments set their setup with `-p` arrays: ISMIP7 Antarctica
  `ismip7_ant` (`scripts/ismip7/params.sh`), Greenland K24 `k24_params`
  (`scripts/ismip7/opt_grl_k24.sh`), 1pctCO2 `pct_grl`, `pct_ant`
  (`scripts/1pctCO2/params.sh`), TIPMIP `tipmip_grl` (`scripts/tipmip/params.sh`),
  nudge `nudge_ant` (`scripts/nudge/run_nudge_antarctica.sh`).
- **Forcing data:** `input/esm/esm_*.nml` (`[esm] par_file`).
- **Scripts:** `scripts/ismip7/` (spin-ups `opt_*.sh`, workflows `run_ismip7_*.sh`;
  see [ISMIP7 projections](experiment-ismip7.md)), `scripts/1pctCO2/` (see
  [1pctCO2](experiment-1pctco2.md)), `scripts/tipmip/` (see
  [TIPMIP Greenland](experiment-tipmip.md)).

## Configuration

- **Run phase and calendar.** `[ctrl] run_step` names the group holding the
  timeline (`"spinup"` or `"transient"`), with `calendar = True` and
  `calendar_ref = 2000.0`: the timeline is in calendar years. The same group holds
  the ESM periods (`time_ref`, `time_hist`, `time_proj`, `time_esm_ref`) and the
  variability switches (`clim_var`, `clim_seed`).
- **Experiment.** `[esm]`: `par_file` (the ESM data configuration in `input/esm/`),
  `experiment` and `esm_name` (the `{experiment}` and `{gcm}` of its file paths),
  `use_esm` / `use_var` / `use_hist` / `use_proj`, `basins` (the basins of the
  ocean extrapolation), and the physical parameters `lapse` (annual mean and
  summer lapse rate; the monthly lapse rate follows a seasonal cosine), `f_p`
  (precipitation scaled by `exp(f_p dT)` with the temperature change `dT` of the
  elevation correction), `f_ocn`, `f_polar`, `dT_threshold`.
- **Anomalies.** In the spin-up the backend applies the reference climatology
  only. In a transient run (`run_step = "transient"`) with `use_esm = True` and
  `experiment` other than `"ctrl"`, it reads the ESM groups of the `par_file` and
  adds anomalies relative to the mean over `time_esm_ref`: up to `time_hist(2)`
  from the historical groups, from `time_proj(1)` on from the projection groups.
  Air temperature and smb anomalies are differences, precipitation a ratio, the
  ocean anomalies are annual means taken at the shelf base. 1pctCO2 and TIPMIP
  point the reference and projection groups at the same file. With `use_esm =
  False`, `experiment = "1pctCO2"` instead applies a spatially uniform anomaly (a
  linear warming of 17 K from 1995 to 2300, capped at `dT_threshold`, times
  `f_polar`; the ocean `f_ocn` times that); `f_ocn`, `f_polar` and
  `dT_threshold` act only there.
- **Forcing on other grids.** A field group of the `par_file` is normally read on
  the climate grid. With `remap` (a kernel: `bilinear`, `shepard`, `nn`, `quadrant`,
  `con`) and `grid_src` (a name for the file's grid, e.g. `"{gcm}-atm"`) it is read
  on its own grid and remapped online: a projected grid with a cdo description
  `grid_<grid_src>.txt` next to the file (as FesmData writes it) or in `maps/`, a
  regular lon-lat grid from its 1D axes (the map is cached in `maps/`), or a
  curvilinear grid (e.g. a tripolar ocean) from the 2D lon/lat named in the
  variable's `coordinates` attribute. `fill = "nearest"` fills the missing
  values from the nearest valid cell. The 1pctCO2 configs read the global CMIP
  fields this way, and the Greenland configs the ORAS5 reference ocean, which
  FesmData provides on GRL-16KM and GRL-32KM only, from GRL-16KM.
- **Surface mass balance.** `[comps] surface_method = "climate"` takes the ESM's
  own SMB (reference + anomaly, corrected from the present-day surface with the SMB
  elevation gradient); `"chion"` (the default) or `"smbpal"` computes it from the
  ESM temperature and precipitation.
- **Update cadence.** `[comps] dt_clim = 1.0`: the forcing is updated every step.
  The spin-up scripts pass `comps.dt_clim=10`: the forcing of a spin-up is a fixed
  climatology, so only the geometry correction changes between updates.
- **Cold start.** `[sim] init_kill_shelves` (no ice where the present-day bed is
  ocean) and `init_time_thrm` (equilibration with topography fixed), shared with
  the other configurations.

## What the backend supplies

- the atmosphere (`tas`, `pr`) for chion or smbpal, or the surface mass balance directly;
- the ocean **at the shelf base** (`T_shlf`, `S_shlf` and their anomalies),
  interpolated from the reference ocean with the marine-shelf parameters, so
  `step_shelf` passes it straight to the marine shelf;
- subglacial discharge `Qd`, landed on Yelmo by `couple_to_yelmo`.

## Output

The shared per-module files (`yelmo.nc`, `isos.nc`, `mshlf.nc`, `chion.nc`, ...),
plus, with `[output] write_clim`, `esm.nc` (the ESM fields on the climate grid) and
`esm_ts.nc` (the forcing averaged over the ice and the floating ice), and, with
`[output] write_cmip` (every `dt_cmip`), the CMIP-formatted `yelmo_cmip.nc` and
`yelmo_ts_cmip.nc`.
