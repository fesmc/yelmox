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
- **Config:** `yelmox/yelmox_esm_*.nml` (Antarctica ISMIP7, Greenland, 1pctCO2, TIPMIP).
- **Scripts:** `scripts/ismip7/` (spin-ups `opt_*.sh`, workflows `run_ismip7_*.sh`;
  see [ISMIP7 projections](experiment-ismip7.md)), `scripts/tipmip/`, `scripts/1pctCO2/`.

## Configuration

- **Run phase and calendar.** `[ctrl] run_step` names the group holding the
  timeline (`"spinup"` or `"transient"`), with `calendar = True` and
  `calendar_ref = 2000.0`: the timeline is in calendar years. The same group holds
  the ESM periods (`time_ref`, `time_hist`, `time_proj`, `time_esm_ref`) and the
  variability switches (`clim_var`, `clim_seed`).
- **Experiment.** `[esm]`: `par_file` (the ESM data configuration in `input/esm/`),
  `experiment`, `esm_name`, `use_esm` / `use_var` / `use_hist` / `use_proj`, and
  the physical parameters `lapse`, `f_p`, `f_ocn`, `f_polar`, `dT_threshold`.
- **Forcing on other grids.** A field group of the `par_file` is normally read on
  the climate grid. With `remap` (a kernel: `bilinear`, `shepard`, `nn`, `quadrant`,
  `con`) and `grid_src` (a name for the file's grid, e.g. `"{gcm}-atm"`) it is read
  on its own grid and remapped online: a regular lon-lat grid from its 1D axes (the
  map is cached in `maps/`), or a curvilinear grid (e.g. a tripolar ocean) from the
  2D lon/lat named in the variable's `coordinates` attribute. The 1pctCO2 configs
  read the global CMIP fields this way.
- **Surface mass balance.** `[comps] surface_method = "climate"` takes the ESM's
  own SMB (reference + anomaly, corrected from the present-day surface with the SMB
  elevation gradient); `"smbpal"` computes it from the ESM temperature and
  precipitation.
- **Update cadence.** `[comps] dt_clim = 1.0`: the forcing is updated every step.
- **Cold start.** `[sim] init_kill_shelves` (no ice where the present-day bed is
  ocean) and `init_time_thrm` (equilibration with topography fixed), shared with
  the other configurations.

## What the backend supplies

- the atmosphere (`tas`, `pr`) for smbpal, or the surface mass balance directly;
- the ocean **at the shelf base** (`T_shlf`, `S_shlf` and their anomalies),
  interpolated from the reference ocean with the marine-shelf parameters, so
  `step_shelf` passes it straight to the marine shelf;
- subglacial discharge `Qd`, landed on Yelmo by `couple_to_yelmo`.

## Output

The shared per-module files (`yelmo.nc`, `isos.nc`, `mshlf.nc`, `smbpal.nc`, ...),
plus, with `[output] write_clim`, `esm.nc` (the ESM fields on the climate grid) and
`esm_ts.nc` (the forcing averaged over the ice and the floating ice), and, with
`[output] write_cmip` (every `dt_cmip`), the CMIP-formatted `yelmo_cmip.nc` and
`yelmo_ts_cmip.nc`.
