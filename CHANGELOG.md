# Changelog

All notable changes to YelmoX are recorded here. Each version corresponds to an
annotated git tag. Dates are release (tag) dates.

## [Unreleased]

### Changed
- Follows yelmo dev with the renumbered `ytherm.qb_method` (1: faces, 2: faces to
  quadrature nodes, 3: simple stagger, 4: quadrature). The par files keep
  `qb_method = 2`, which now selects the energy-consistent "faces to quadrature
  nodes" method (was quadrature), so results change. Par files drop
  `ytopo.surf_gl_method`, `ytopo.margin2nd` and `ydyn.ssa_beta_max` (removed in
  yelmo); `input/` yelmo copies re-synced. Requires that yelmo dev.
- `input/`: yelmo input copies re-synced with yelmo kryos-init (`yelmo_defaults.nml`:
  `yelmo.mask_border`; `yelmo-variables-ydyn.md`: `H_ice_solv`, `f_ice_solv`).
  Requires a yelmo with these keys.
- **New `[domain]` group defines the domain** (`[domain_north]`/`[domain_south]`
  in bipolar): `name`, the grid of every component (`grid_hub`, `grid_ice`,
  `grid_isos`, `grid_clim`, `grid_smb`, `grid_mshlf`; blank = default) and the
  hub's topography and code masks (`topo_path`, `topo_names`, `regions_path`/`_var`,
  `basins_path`/`_var`; a blank mask path gives a mask of 1). It replaces
  `[htopo]`, the `grid_*` keys of `[coupling]` and `domain`/`grid_name` in
  `[yelmo]` (Yelmo gets them from `[domain]`). All par files and run scripts are
  migrated (`yelmo.grid_name`/`htopo.grid_name` -> `domain.grid_ice`/`grid_hub`,
  `coupling.grid_*` -> `domain.grid_*`); results are unchanged. The hub grid
  no longer tracks the Yelmo grid: `grid_ice = ""` tracks the hub instead.
- The hub also reads the bed roughness `z_bed_sd` (4th `[domain] topo_names`
  entry; `""` = none, 0). Par files take the name from `yelmo_init_topo`, except
  where that was not a bed-roughness field (`bed`, `bed_bedmap3`, `H_ice`,
  `none`), which become `""`. SRG keeps `z_bed_err`.
- **Yelmo is populated from the domain.** Its grid comes from
  `maps/grid_<grid_ice>.txt`, and the hub topography, remapped conservatively to
  `grid_ice`, is both its initial topography and its present-day reference
  (`H_ice_ref`, `z_bed_ref`, optimization target). Par files drop `yelmo.grid_path`,
  `yelmo_init_topo.init_topo_load/path/names` and `yelmo_data.pd_topo_load/path/names`.
  Unchanged where `grid_ice = grid_hub`; where they differ (yelmox Antarctica,
  bipolar south: Yelmo 32 km, hub 16 km) Yelmo starts from the remapped 16 km
  topography instead of its own 32 km files. `f_grnd_pin` (diagnostic) changes
  where `z_bed_sd` was read from a non-roughness field. Requires yelmo kryos-init.
- The hub fills the gaps of its topography file (missing values, e.g. outside
  the coverage of the ISMIP7 obs files): `H_ice = 0`, `z_bed` from the nearest
  valid cell (fesm-utils `fill_nearest`), `z_srf` from `z_bed` and `H_ice` at sea
  level 0, and `z_bed_sd = 0`. The counts are logged. Before, the raw fill
  values (-9e33) reached the remaps; Yelmo's own reads set them to -9999.
- **Ice mask and named regions from the domain.** `[domain]` gains
  `sectors_path`/`_var` (a third code mask), `ice_codes_mode` (`all`, `include`,
  `exclude`) with `ice_codes` (codes of `regions`: where ice is allowed), and
  `region_names`/`region_mask`/`region_codes` (named regions for 1D output,
  `yelmo_ts_<name>.nc`). The hub's `regions` and `basins` (nearest neighbour to
  `grid_ice`) are now Yelmo's too, one set for every component, and the ice mask
  is passed to `yelmo_init`; `[yelmo_masks]` keeps no keys. This replaces
  `libs/ice_sub_regions.f90` and the per-domain masks in `domain_regions_init`
  and in Yelmo. Par files: Antarctica `exclude 2.0`, APIS/WAIS/EAIS = sectors
  3/1/2 of `BASINS-nasa mask_regions`; Greenland `include 1.3 1.11 1.0`;
  Laurentide `exclude 1.30`, Hudson = regions 1.12, `yelmo.mask_border = "none"`;
  North `exclude 1.0`; others `all`. Physics is
  unchanged where `grid_ice = grid_hub`; Yelmo's `basins`/`regions` output
  changes where its own files differed (Antarctica `basin_reese` -> `basin`).
- **Relaxation to the reference from the domain.** `[domain]` gains
  `relax_codes_mode` (`none`, `all`, `include`, `exclude`), `relax_codes` (codes of
  `regions`) and `relax_tau`: Yelmo's `tau_relax` is `relax_tau` where selected and
  -1 (free) elsewhere, used with `ytopo.topo_rel = -1`. It replaces the Patagonia
  case of `domain_regions_init`, which no config reached since the domain was
  renamed SRG. `yelmox_SRG.nml`: `exclude 1.0`, 50 yr (the icefield evolves freely,
  the rest relaxes, as the Patagonia case did); all other par files `none`.
- The domain type `ice_domain` is renamed `kryos_domain`, in line with the
  Kryos naming of the cryosphere-component framework.
- `libs/yelmox_domain.f90` is split, by concept, into `kryos` (domain type,
  configuration, init, `remap`), `kryos_regions` (region-specific masks and
  physics), `kryos_coupling` (`step_*`, `couple_*_to_yelmo`), `kryos_startup`
  (cold start, restart bundles), `kryos_output` and `kryos_forcing` (`tsforcing`).
  Code is moved unchanged; drivers import each name explicitly.
- Renamed coupling primitives: `step_optimize` -> `step_spinup_tuning` (it also
  ramps the relaxation timescale), `refresh_htopo` -> `refresh_hub`,
  `domain_update_smb` -> `step_smb`.
- `couple_to_yelmo` assembles the Yelmo boundary state as its own step, called
  by the drivers before `step_icesheet` (which no longer runs the couplers).
  The Greenland NEGIS friction update now sees the bedrock of the current step
  (`use_negis = True` only; no config sets it).
- `step_climate` no longer runs the surface mass balance; drivers call
  `step_smb` right after it.
- `yelmox` and `yelmox_bipolar` write the per-step coupling sequence out in the
  time loop; `yelmox_step` and the bipolar `advance_isostasy`/`advance_dynamics`
  wrappers are gone.
- `step_climate` and `domain_startup` take the transient forcing object (`tsf`)
  as one optional argument instead of `dTa`/`dTo`/`dSo`; `update_climate` applies
  its anomalies only when it is active.
- `domain_ctl` grid names: `grid_name` -> `grid_hub` (the hi-res hub),
  `grid_yelmo` -> `grid_ice` (Yelmo).
- Cold starts made consistent across drivers. `yelmox_esm` and `yelmox_rembo`
  now set up isostasy through the shared `domain_init_isostasy` (conservative
  ice-load coarsening + isostasy reference check). The optimisation's
  `cb_ref = cf_init` guess is set before `yelmo_init_state` in every driver
  (was after it in `yelmox`/`yelmox_bipolar`), so results of `opt` cold starts
  change slightly. `opt.cf_init` must be > 0; the "negative uses cb_tgt" rule
  (only implemented in `yelmox_rembo`) is gone.
- `yelmox_rembo`: `greenland_init_marine_H` applies the shared rule (H = 800 m
  where H < 600 m and z_bed > -500 m) instead of H×1.2, and the driver no longer
  shrinks `dtt`/`dtime_emb` during a tsgen ramp.
- `yelmox_bipolar`: the OBM restart (`obm_restart.nc`) is written into the
  run-root restart bundle and read back from `[ctrl] restart_bsl`;
  `&nautilus use_restart/restart` removed.

- Driver time loops (`yelmox`, `yelmox_bipolar`, `yelmox_esm`, `yelmox_rembo`):
  output and restarts are written at the top of the loop for the current time
  (`time_init` on the first pass), then the loop exits once finished, else
  `tstep_update` advances the time and the domain is stepped. Each output call
  appears once; the final state and a final restart bundle are written on the
  last pass, also when the tsgen kill switch trips. Requires fesm-utils with
  `tstep_update` advancing on every call.
- The zero-length first step at `time_init` is gone. It also ran an extra
  basal-friction optimisation update (and, in `yelmox_bipolar`, an extra OBM
  step) at every start, so results of cold starts change slightly.
- Restart cadence: `[coupling] dt_restart` is replaced by a `[tm_rst]` timeout
  group (`method = none|const|file|times`, `dt`, `file`, `times`); a restart
  bundle is always written at `time_end`. Par files updated (`dt_restart = X`
  -> `method = "const"`, `dt = X`). In `yelmox_bipolar` the cadence is shared
  by both domains.
- Timeline set up with fesm-utils `tstep_init(ts, path_par, group, dtt, ...)`
  (namelist form); `tstep_due` renamed `cadence_due` (only `dt_clim` and the
  esm CMIP output use it).

- `input/`: yelmo input copies re-synced with yelmo:dev (`yelmo_defaults.nml`,
  `yelmo-variables-{ydyn,ytopo}.md`; added `elsa_defaults.nml`,
  `tracer_defaults.nml`). Par files drop `ydyn.scale_T` / `ydyn.T_frz`
  (no longer read by yelmo).
- `&opt`: `use_yelmo_cf_min`, `opt_cf_min` and `cf_min` removed from all par
  files (removed from yelmo's optimizer; the `cf_ref` floor is `ytill.cf_min`).
- `scripts/1pctCO2/opt_ant.sh`, `scripts/ismip7/opt_{ant,grl}.sh` ported from
  `esm-legacy` to `esm` (par files `yelmox_esm_Antarctica_1pctCO2.nml`,
  `yelmox_esm_Antarctica_ismip7.nml`, `yelmox_esm_Greenland.nml`; keys mapped:
  `coupling.equil_method`, `ytill.cf_min`, `ycalv.tau_ice_flt`,
  `yhyd.bkt_N_closure`/`marine_p`).

### Fixed
- `yelmox_esm_Antarctica.nml`, `yelmox_esm_Antarctica_nudge.nml`: `&ghf` lacked
  `obs_err_name` and `f_stdev` (startup stopped on the nml read).
- `yelmox_esm_Antarctica.nml`: topography and geothermal heat flux read from
  the current `{grid}_TOPO-BedMachine.nc` and `{grid}_GHF-HR24.nc` (the old
  `TOPO_BedMachineAntarctica-v3` and `GHF-M17` files are no longer in
  `ice_data`).
- `input/esm/esm_ant_1pctCO2.nml`: meltMIP OI ocean reference read from
  `ice_data/ISMIP7/` (shared with ISMIP7), not a separate `1pctCO2/` copy.
- `scripts/1pctCO2/opt_ant.sh`: `resolution` now sets `yelmo.grid_name`
  (it only named the output folder).
- `yelmox`, `yelmox_bipolar`: after restoring a restart bundle the climate/smb
  and marine-shelf forcing are rebuilt before the first step (as in
  `yelmox_esm`/`yelmox_rembo`); the bundle does not hold them.
- `esm_forcing`: transient ocean anomalies (`dto`, `dso`) were NaN/garbage for
  depth-less ocean forcing (Greenland ISMIP7 `tf`/`so`, 2D monthly): the months
  were averaged over axis 4, which holds them only for 3D fields (out-of-bounds
  reads; the transient stopped at its first step). The annual mean now uses
  fesm-utils `varslice_sub_mean`; the hist/proj `to`/`so` blocks share
  `esm_ocean_anomaly`. Annual 3D ocean (Antarctica) unchanged. Requires
  fesm-utils with `varslice_sub_mean`.
- `esm_forcing_init`: a transient ocean field whose layout (rank, extent, depth
  levels) differs from its ESM reference stops with an error.
- `yelmox_esm` Greenland: `dT_shlf` was `T_shlf + dto`, i.e. the absolute shelf
  temperature (K) plus `dto` again, so with `tf_method=2` (`bmb_method="anom"`)
  `tf_shlf` was ~274 and shelf melt hundreds of m/yr. `dT_shlf`/`dS_shlf` are
  now the anomalies `dto + dto_var` / `dso + dso_var` for every domain, and the
  Greenland-only override of `tf_method` is gone (the par files set it).
  Antarctica (`tf_method=1`) is unchanged. Greenland ESM spin-ups need re-running.
- `yelmox_esm_Antarctica_ismip7.nml`: `&spinup time_ref` was 1961-1990, outside
  the 1979-2022 axis of the RACMO2.3 reference climatology, so the reference
  `t2m`/`pr` were missing values (-9999) and the ice sheet melted away in the
  first step. Now 1985-2014, as in `&transient`. ISMIP7 Antarctica spin-ups
  need re-running.
- `yelmox_esm_Greenland{,_1pctCO2,_tipmip}.nml`: `&itm` had the Antarctic ITM
  parameters (`itm_c=-55`, `itm_b=3`, `itm_lat0=-60`, `alb_ice=0.70`), so the
  latitude-adjusted `itm_c` was about +340 W m-2 over Greenland and the SMB was
  strongly negative everywhere. Now the Greenland values of
  `yelmox_Greenland.nml` (`-45`, `-2`, `65`, `0.4`). Greenland ESM spin-ups need
  re-running.
- `yelmox_esm` output: `pr_ann` (2D) and `pr_1d` (time series) were the
  precipitation in mm/d times 1e-3, labelled m/a, with only the January anomaly
  factor. Now the annual mean of `pr*dpr`, in mm/d.
- `yelmox_esm` output: the monthly anomalies entered `t2m_ann`/`t2m_sum` (small
  2D file) and `t2m_1d`, `dt_1d`, `dpr_1d`, `dt_var_1d`, `dpr_var_1d` with January
  only, and `t2m_sum` (2D file) with the DJF mean (scaled by 0.333) also in the
  north. Now annual means, and summer means of the hemisphere (DJF south, JJA
  north; new `esm_summer_mean`, also used for `esm%t2m_sum`).
- `smbpal`: the state fields are zeroed at allocation. With `abl_method="pdd"`
  the albedo `alb_s` was never set, so `smbpal_restart.nc` held uninitialized
  memory (differing between builds). Results are unchanged.

### Removed
- `libs/simpleclim.f90` (empty stub) and the unreachable `"const"` branch of
  the LGM-north cold start.
- `timeline_init` (replaced by `tstep_init`) and `domain_ctl%dt_restart`.
- Legacy single-grid programs (`<flavor>/legacy/`, `make <flavor>-legacy`) and
  retired flavors (`retired/`: `yelmox_ismip6`, `yelmox_nahosmip`,
  `yelmox_rtip`), with their par files, make targets, runme aliases and
  `scripts/ismip6-2300.md`. They no longer ran against yelmo:dev.

## [v2.3] - 2026-07-15

### Added
- Multigrid driver framework: all flavors (`yelmox`, `yelmox_bipolar`,
  `yelmox_esm`, `yelmox_rembo`) rebuilt on the shared `ice_domain` type and
  `step_*` coupling primitives in `libs/yelmox_domain.f90`; per-flavor
  documentation pages.
- Transient time-series forcing (`tsgen`/`tsforcing_class`): single forcing value
  mapped onto per-channel anomalies, forcing-increment (`Δf`) restarts, kill
  switch, and 1D diagnostics; wired into the snapclim, ESM, and REMBO drivers.
- `snapesm` climate backend: backend-agnostic `dom%clim` adapter
  (`yelmox_climate`) so an ESM/varslice climate can replace snapclim.
- `yelmox_esm`: ISMIP7 spin-ups (Greenland + Antarctica), TIPMIP Greenland
  stabilisations, and 1pctCO2 forcing-only scaffolds; annual-mean transient ocean
  via `varslice_nsub`.
- Namelist database in the compact fesm-utils:dev `varslice` format, with a
  database-namelists documentation page.

### Changed
- `yelmox_esm`: `esm` module/source renamed to `esm_forcing` (avoids ifx name
  clash); forcing loaded through shared `varslice_init_nml` with
  `{gcm}`/`{experiment}` substitutions.
- Sync with yelmo:dev: `ytrc` tracer subsystem, `var_io` tables, and default
  parameters.
- `yelmox_rembo`: driver-contained routines take explicit arguments instead of
  relying on host association, matching `yelmox_esm`/`yelmox_bipolar`.
- Output layout: drop grid suffix, add `yelmo_sm.nc`, htopo mask-load flags.

### Removed
- Retired the `yelmox_ismip6`, `yelmox_nahosmip`, and `yelmox_rtip` flavors to
  `retired/`. They still compile (`make yelmox_ismip6` / `yelmox_nahosmip` /
  `yelmox_rtip`, each printing a retirement notice); prefer `yelmox` or
  `yelmox_esm` for new work.
- Dropped pre-configme / pre-runme tooling.

### Fixed
- `var_io` tables (`input/yelmo-variables-{ydata,ytrc}.md`) synced to the
  yelmo:dev isochrone dimension rename (`age_iso`→`time_iso`,
  `pd_age_iso`→`pd_time_iso`). The stale copies crashed restart writing
  (`nf90_inq_dimid`), which affected every run.

## [v2.2.2] - 2026-06-24
- yelmox: added bsl ts writing.

## [v2.2.1] - 2026-06-24
- Config and small bug fixes.

## [v2.2] - 2026-06-18
- Yelmo default parameters added. `yhyd` section (with bug fix to bucket units).
  Added yelmo-config tool for inspecting and comparing parameters. yelmox: new
  capabilities to support ISMIP7 activities.

## [v2.1.3] - 2026-06-15
- Added support for ISMIP7 simulations GrIS+AIS.

## [v2.1.2] - 2026-06-15
- Added `with_isostasy` parameter for yelmox - now complete.

## [v2.1.1] - 2026-06-15
- Added `with_isostasy` parameter for yelmox.

## [v2.1] - 2026-06-11
- Conversion to FastHydrology for yelmo basal hydrology.

## [v2.0.6] - 2026-06-02
- Added missing yelmo variables in io tables.

## [v2.0.5] - 2026-06-01
- Added missing isos parameter, removed obsolete.

## [v2.0.4] - 2026-05-30
- yelmox/yelmo work with configme v0.6.6+ and runme v0.5.8+.

## [v2.0.3] - 2026-05-30
- yelmox/yelmo work with configme v0.6.2+ and runme v0.5.8+.

## [v2.0.2] - 2026-05-24
- Tagged version consistent with yelmo:v2.0.2, with central runme and configme setup.
