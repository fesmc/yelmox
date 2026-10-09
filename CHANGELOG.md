# Changelog

All notable changes to YelmoX are recorded here. Each version corresponds to an
annotated git tag. Dates are release (tag) dates.

## [Unreleased]

### Added
- chion BESSI in yelmox (`[chion] model = "bessi"`): the host supplies its surface shortwave
  (TOA insolation times `[surface_chion] trans_sw`), constant `wind_speed` and `rel_hum`, the
  air pressure from the surface elevation (barometric, annual-mean air temperature) and the
  solar longitude; longwave from chion's own parameterization. ITM results unchanged.
- 1pctCO2 (Greenland, Antarctica) and TIPMIP (Greenland) finished: ice sheet and isostasy
  on, ISMIP7-like optimization spin-up; TIPMIP stabilisations branch off the ramp's
  restarts (years 109, 232). 1pctCO2 reads the global CMIP fields of `[esm] esm_name`
  (default MPI-ESM1-2-LR) with online remapping (fesm-utils `varslice` `remap`/`grid_src`).
- `[ghf] convert_ghf_units` (W/m2 -> mW/m2); missing (negative) GHF values are filled
  from the nearest valid value.
- Greenland ISMIP7 setup with the K24 basal hydrology: `yelmox/yelmox_esm_Greenland_k24.nml`
  (`yelmox_esm_Greenland.nml` with `method_transport = 1`, so K24 sets N_eff, and
  `k24_sliding_law = 4`, so K24 takes Yelmo's basal stress and its water source includes
  the frictional heat that `q_T` already contains) and the spin-up `scripts/ismip7/opt_grl_k24.sh`.
- Antarctic paleo setup (32 km): `yelmox/yelmox_Antarctica_paleo_spinup.nml`
  (15 kyr optimization spin-up) and `yelmox_Antarctica_paleo_lgp.nml` (-130 kyr to
  +10 kyr, climate from the glacial index `input/alpha_combined_125kyr_interp.dat`,
  sea level from `sealevel_rohling_450kyr.dat`, ages traced with elsa), with run
  scripts in `scripts/ant-paleo/`. Ported from the old single-grid par files; the
  transient now runs on relative time so that
  sea level follows the record (it stayed at present day before).
- snapesm writes its state to `snap.nc` (`[output] write_clim`): the driving
  indices `idx_<name>`, `z_srf`, monthly `tas`/`tsl`/`pr`, `ta_sum`, `tsl_ann`
  and the ocean profiles `to_ann`/`so_ann` (before only `t2m_ann`, `pr_ann`).
- snapesm restarts: restart bundles hold the state of each driving index
  (`snapesm_idx_<name>_restart.nc`, with the snapesm configuration as
  attributes), restored on restart; bundles without them keep the cold-start
  indices.
- `scripts/ismip7/`: the ISMIP7 optimization spin-ups (`opt_ant.sh`, `opt_grl.sh`,
  and L. Gutierrez Gonzalez's `opt_grl_ismip.sh`) for the current `yelmox_esm`.

### Changed
- chion is the default surface model: `[comps] surface_method = "chion"` in every ITM par
  file (Greenland, Antarctica paleo, SRG, pd_Greenland, all ESM except `nudge`, bipolar
  north) and the kryos built-in default. Each gets `[surface_chion]` + `[chion]` with the
  host values of its own `[smbpal]`; `[smbpal]` stays, so `surface_method = "smbpal"`
  switches back. The PDD par files (Antarctica, North, LIS, Pyrenees, pd_Antarctica,
  esm_Antarctica_nudge, rembo, bipolar south) keep smbpal.
- Par files: `ydyn.ssa_lis_opt_energy` CG with SSOR preconditioner and `-tol 1.0e-3` (was
  Jacobi, `1.0e-4`, which hit `-maxiter 200` in almost every ANT-8KM solve; the result was
  ~1e-3 converged anyway). ANT-8KM: 1/179 solves at the limit, 52 instead of 200 iterations.
- smbpal ITM runs point by point with OpenMP; the daily forcing and insolation are prepared
  once per call (once for the whole 100-yr cold-start snowpack equilibration). ANT-8KM
  ISMIP7 10-yr test: 22.5 -> 8.0 min (16 threads). Equal to the serial version to
  single-precision round-off. Daily smbpal output (`file_out_day`) removed.
- `scripts/check_par_nml.py` follows yelmo's defaults file: it reports keys of yelmo groups
  missing from `input/yelmo_defaults.nml` (what `nml_validate` stops on) and duplicate
  groups/keys, no longer omitted parameters; `--dead` also reports non-yelmo keys that no
  source reads. Checks `yelmox/*.nml` and `yelmox_bipolar/*.nml` by default.
- **ISMIP7 geothermal heat flux**: the ObsISMIP7 GHF (W/m2) was never converted and,
  with `obs_err_name` = the GHF itself, clamped to 0.1 mW/m2 everywhere. Fixed
  (`convert_ghf_units = True`, `obs_err_name = "none"`); ISMIP7 spin-ups must be redone.
- `[esm]` subglacial discharge is read only if the sgd group's `filename` is not
  `"none"` (Qd = 0), instead of always for Greenland; `[esm] grid_src` removed (unused).
- `climate_init`/`esm_forcing_init` take the climate `grid_class` instead of its name and size.
- `kryos` builds no Yelmo<->hub maps when the grids are the same.
- `sed.method` is read (was uninitialised); `use_obs` keys replaced by `method`.
- Greenland par files use `TOPO-M17-v5`; pd_Greenland/pd_Antarctica rebuilt from
  yelmox_Greenland/Antarctica (present-day climate, 50 kyr); production run lengths
  (Antarctica 15 kyr, ESM spin-ups 15/20 kyr) instead of bring-up values.
- Par files: dead keys removed (timeline `use_hist`/`use_proj`, `snap_pd.is_ref`,
  `opt.tf_time`, `opt.cf_ref_wais`, `marine_shelf.basin_name/basin_bmb`,
  `isos.dt_prognostics`, Pyrenees `*_L21`); `opt.basin_fill` added where missing;
  personal paths removed. `input/snap_index` (link into a home directory) removed.
- `.runme/info.json`: the REMBO namelist is staged for `yelmox` runs (its key was
  still the removed `yelmox_rembo` executable).
- Build: yelmox objects depend on `libfesmutils.a`, `libyelmo.a` and `libisostasy.a`
  (`climate_rembo.o` also on `librembo.a`), so a type change in a dependency rebuilds
  them instead of leaving a stale layout (fesmc/FastHydrology#10). Needs fesm-utils,
  yelmo, FastIsostasy and rembo1 dev with the archive as a file target, else every
  `make` rebuilds yelmox. FastIsostasy and rembo1 objects depend on `libfesmutils.a`.
- `ytopo.gl_sep = 3` instead of 2 in `yelmox_esm_Antarctica*.nml` and the ISMIP7/1pctCO2
  optimization scripts (`opt_ant.sh`, `opt_grl*.sh`): yelmo (b7049784) removes `gl_sep = 2`,
  which gave `f_grnd = 0` to cells grounded at their centre next to deep ocean (GRL-16KM
  spin-up killed at t = 1 yr). Needs fesm-utils cc3f719.
- Follows yelmo dev (f78eaa91): `input/yelmo_defaults.nml` gains `yelmo.pc_rho_max` (2)
  and `ytherm.gl_temperate` (True); `ydyn.slide_T`/`gamma_T`/`lambda_min` are replaced by
  `frz_scale`/`frz_efold`/`frz_min` (3 K, 1e-3; sliding-speed e-fold, was ~0.2 K for
  `beta_q = 0.2`) and `ytherm.use_strain_sia` by `strain_heating = "full"`. All par files
  take `pc_tol = 1` (was 5), `bkt_floating_mode = 0` (was 1) and the default `de_max = 100`
  (was 0.5). Ice-sheet results change.
- Follows yelmo dev (6fad9a11): `&yhyd` K24 keys renamed as in FastHydrology dev
  (fesmc/FastHydrology#14), e.g. `k24_ub_hook` -> `k24_N_ub_coupled`, `k24_eta_w` ->
  `k24_water_viscosity` (full list in the yelmo changelog); new `k24_kappa_z_hard` /
  `k24_kappa_z_soft`. Values unchanged. Needs FastHydrology dev 8e681d0 or later.
- `[opt] opt_cf` and `opt_tf` are methods instead of switches: `opt_cf = "none" | "L21"`,
  `opt_tf = "none" | "L21" | "L21-points"`. `"L21"` optimizes one `tf_corr` per basin
  (`optimize_tf_corr_basin`, `tf_basins`), `"L21-points"` each point (`optimize_tf_corr`,
  `tf_sigma`, `basin_fill`; the method used so far). Par files: `True` -> `"L21"`
  (`opt_cf`) / `"L21-points"` (`opt_tf`), `False` -> `"none"`. Logical values stop the
  model. Needs yelmo with string `opt_cf`/`opt_tf` (`libs/ice_optimization.f90`).
- The topography relaxation of a spin-up is its own switch, `[sim] relax`, with
  the group `[relax]` (`topo_rel`, `tau1`, `tau2`, `time1`, `time2`, `m`; was
  `[opt] rel_tau1/2`, `rel_time1/2`, `rel_m` with `topo_rel = 4` fixed), applied
  by `step_relax` before `step_optimize`. After `time2` the `[ytopo]` values of
  `topo_rel` and `topo_rel_tau` apply again (was `topo_rel = 0`). `step_optimize`
  only optimizes `cb_ref` and `tf_corr`, so `opt` can run without relaxation.
  The par files set `relax` as their `opt`; the scripts pass `sim.relax=True`
  with `sim.opt=True`. Needs yelmo with `relax_params` (`libs/ice_optimization.f90`).
- Components named by the boundary they supply: `surface` (mass balance +
  temperature) and `shelf` (shelf-base melt + temperature). Keys:
  `with_marine_shelf` -> `with_shelf`, `smb_method` -> `surface_method` (values
  unchanged), `grid_smb` -> `grid_surface`, `grid_mshlf` -> `grid_shelf`,
  `write_smb` -> `write_surface`, `write_mshlf` -> `write_shelf`; new
  `[comps] with_surface` (the surface was switched by `with_climate` before).
  `with_surface` and `with_shelf` need `with_climate` (else the run stops).
  Routines: `step_smb` -> `step_surface`, `step_marine_shelf` -> `step_shelf`,
  `couple_smb_to_yelmo` -> `couple_surface_to_yelmo`, `couple_marine_to_yelmo`
  -> `couple_shelf_to_yelmo`. Output files keep the model names (`smbpal.nc`,
  `mshlf.nc`).
- `[coupling]` is split into `[comps]` (`with_*`, `climate`, `smb_method`,
  `dt_clim`: which components are active, with which model, how often) and
  `[sim]` (the conditions of the simulation: cold-start ice state, optimization,
  regional modifications); in `yelmox_bipolar` `[comps_<sfx>]`, `[sim_<sfx>]`.
  Renamed in `[sim]`: `equil_method = "none"/"opt"` -> `opt = False/True` (the
  optimization runs within the `[opt]` time windows, on a cold start or a
  restart), `kill_shelves` -> `init_kill_shelves`, `time_equil_thrm` ->
  `init_time_thrm`. Old par files stop with "parameter not found".
- The restart bundle to start from is `[ctrl] restart` (was `[coupling] restart`;
  in `yelmox_bipolar`, `[ctrl] restart_bsl` and `[coupling_<sfx>] restart`). The
  driver reads it and passes it to `domain_startup(dom, ts, bsl, restart, ...)`.
  `yelmox_bipolar` writes one bundle per run: `restart-<kyr>-kyr/` holds the
  shared bsl and obm restarts and each domain in a subfolder named after it
  (was `<domain>/restart-<kyr>-kyr/`).
- Follows yelmo dev (2c3d3449; needs yelmo dev at or after 9e93696d): `input/yelmo_defaults.nml`
  gains `ydyn.ssa_vel_lim_method` (default `"drag"`, a smooth speed-limit drag) and
  `ssa_vel_lim_tau`; all par files take `ssa_vel_max = 10000` (was 5000) and
  `pc_eps = 0.02` (was 1.0). Ice-sheet results change.
- The hemisphere of a domain (seasons, lapse rates in snapclim, snapesm and the
  esm forcing) follows from the latitude of `grid_clim` (south when its mean is
  below 0), no longer from the domain name `Antarctica`.
- `yelmox_rembo` is removed: REMBO runs use `yelmox` with `climate = "rembo"`
  and `smb_method = "climate"`, built with `make yelmox rembo=1` (without it, a
  stub stops the run). REMBO supplies the atmosphere and smb, snapclim the ocean.
  The par files move to `yelmox/yelmox_rembo_Greenland.nml` and
  `yelmox/rembo_Greenland.nml` (`ctrl.time_equil` becomes
  `coupling.time_equil_thrm`, `ctrl.write_restart` is dropped), the run script to
  `scripts/rembo/` (runme `-e yelmox`; the `rembo` alias is gone). Output with
  `[output] write_clim`: `rembo.nc` (`t2m_ann`, `t2m_sum`, `pr_ann` now in mm/a,
  `smb_ann`) and `rembo_ts.nc` (`dT_sum`, the applied `f_now·f_ta`, was `dT_jja`
  = `f_now`; `dT_ann`, `dT_ocn`, `smb_mean`, `aar`; `V_dT` is dropped); the tsgen
  forcing is in `yelmo_ts.nc`.
- With ESM forcing the marine shelf takes the ice-shelf base as `z_srf - H_ice`
  (Yelmo's definition), like the other climates, instead of reconstructing it from
  flotation. ESM results change under the ice shelves.
- `yelmox_esm` is removed: ESM runs use `yelmox` with `climate = "esm"`. Its par
  files move to `yelmox/yelmox_esm_*.nml` (without `[esm] use_smb`), its run
  scripts to `scripts/ismip7/`, `scripts/tipmip/` and `scripts/1pctCO2/` (runme
  `-e yelmox`; the `esm` alias is gone), and scripts set
  `coupling.smb_method = climate | smbpal` instead of `esm.use_smb`.
- Output: the climate backend writes its own file, `[output] write_clim`
  (was `write_snap`): `snap.nc` (snapclim, snapesm) or `esm.nc` (the esm fields:
  temperature, precipitation or SMB anomalies, shelf anomalies, discharge) plus
  `esm_ts.nc` (the forcing means over ice and floating ice, before in
  `yelmo_ts_esm.nc`, now on the climate grid). CMIP/ISMIP-formatted output is a
  general option, `[output] write_cmip` / `dt_cmip` (was `[esm] write_formatted` /
  `dt_formatted`), with its writers in `libs/cmip_output.f90`. Scripts set
  `output.write_cmip`.
- ESM forcing is a climate backend: `[coupling] climate = "esm"` runs
  `esm_forcing` (unchanged) inside `yelmox`, reading `[esm]` and its periods from
  the `run_step` group. The climate products grow: the ocean at the shelf base
  (`T_shlf`/`S_shlf`, used by `step_marine_shelf` in place of depth profiles), the
  surface mass balance (`smb_method = "climate"`: reference smb + anomaly,
  corrected from the present-day surface with the smb gradient; replaces
  `[esm] use_smb`) and subglacial discharge (`Qd`, landed by `couple_to_yelmo`).
  `climate_update` takes the domain geometry on `grid_clim`, the marine-shelf
  parameters and `dtt`. The ESM par files set `dt_clim = 1` (every step) and gain
  `[tsforcing]`/`[tsgen]` (inactive).
- Run control in `[ctrl]`: `run_step` names the group holding the timeline
  (`"ctrl"` = `[ctrl]` itself; ESM `"spinup"`/`"transient"`), `calendar` /
  `calendar_ref` set calendar years and their reference (ESM `True`/`2000`, before
  hard-coded in `yelmox_esm`). Cold start in `[coupling]`, for every driver:
  `kill_shelves` (no ice where the present-day bed is ocean) and `time_equil_thrm`
  (equilibration with topography fixed after `init_method`), applied in
  `domain_init_ice`. They replace the ESM `[spinup]`/`[transient]` keys
  `kill_shelves` / `time_equil`; the fixed-topography equilibration no longer
  requires `equil_method = "opt"`. Scripts set `coupling.kill_shelves`.
- The climate backend is chosen at runtime: `[coupling] climate = "snapclim" |
  "snapesm"` replaces the `make CLIMATE=` switch (one `libs/yelmox_climate.f90`
  holds both backends; `yelmox_snapesm.x` is gone). The driver's transient forcing
  goes to the backend (`climate_update(..., tsf)`), which applies it its own way.
  The bipolar ocean box model takes its air-temperature anomaly from
  `climate_air_anom` and rembo asks `climate_ocean_const`, instead of reading
  snapclim's internals. Par files set `climate` (`"esm"` in the `yelmox_esm` ones,
  not read until ESM becomes a backend); results do not change.
- The domain owns its physical constants: `[domain] phys_const` (e.g. `"Earth"`,
  a group of `input/yelmo_phys_const.nml`) is loaded once and the same record goes
  to the hub, Yelmo (`yelmo_init` `cnst`), isostasy, the marine shelf and
  `smb_simple`. Before, Yelmo loaded them and the others took Yelmo's copy; the
  hub had its own densities. `yelmo.phys_const` now only selects Yelmo's calendar
  year. All par files set `phys_const = "Earth"`; results do not change.
- The hi-res hub keeps its reference geometry (`z_bed_ref`/`H_ice_ref`/`z_srf_ref`).
  On a hub finer than Yelmo, `couple_yelmo_to_htopo` (was `refresh_hub`) adds
  Yelmo's bed displacement and change in ice thickness to it and recomputes
  grounding (from flotation) and the surface on the hub (`htopo_update`), instead of overwriting the hub with
  Yelmo's fields refined bilinearly. On Yelmo's grid the hub still mirrors
  Yelmo. Multigrid runs (e.g. Antarctica, hub 16 km / Yelmo 32 km) change
  through the marine shelf and the climate's surface elevation.
- Builds use OpenMP by default (`openmp ?= 1` in `config/Makefile`); `make <driver>
  openmp=0` builds serial. Regenerate the Makefile with configme to pick it up.
- `input/yelmo_defaults.nml` re-synced with yelmo dev (`ytrc.elsa_restart`).
- Follows yelmo dev: `input/` yelmo copies re-synced (`yelmo.mask_border`, `"auto"`:
  the domain border as before; the capacity basal BC keys of `ytherm`, now the
  default, so results change; the K24 options of `yhyd`). Par files rename
  `yhyd.k24_long_coupling_water = 5` to `k24_coupling_length_kamb86 = 10` (yelmo's
  rename, twice the old value) and set `yhyd.k24_flux_solver = 3` (the taped solver,
  which FastHydrology's default routing scheme now requires; K24 transport is off
  in all par files, so results do not change). Requires that yelmo dev.
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
- Greenland no longer sets `cb_ref = ytill.cf_ref` at the start when
  `ytill.method = -1`: an external `cb_ref` is up to the user (optimization or
  restart). A Greenland cold start with `method = -1` and no optimization (e.g.
  `yelmox_esm_Greenland*.nml` with `equil_method = "none"`) now starts from
  Yelmo's `cb_ref = 1` fallback instead of 40.
- `scale_glacial_smb` and `use_negis` apply to any domain, not only to one named
  Greenland. `scale_glacial_smb = True` needs a `[glacial_smb]` group (`dt_lgm`,
  `lat_lim`, `fac_lim`; were fixed at -8 K, 55°N, 0.9); `[negis]` gains
  `basin_centre`/`basin_south`/`basin_north` (were fixed at 9.1/9.2/9.3). No par
  file sets either switch.
- **Cold-start ice state from `[coupling]`**, the same in every driver
  (`domain_init_ice`): `init_marine_H` (was `greenland_init_marine_H`), then
  `init_method` = `none`, `equil` (`init_equil_time`), `recon` or `recon_ref`
  (`recon_path`, `recon_var`, `recon_codes`). It replaces the startup chosen by
  domain name (and, for Laurentide/North, by `tstep_method`) and the own
  equilibrations of `yelmox_esm` (1 yr, spin-up only) and `yelmox_rembo` (10 yr).
  Par files keep their behaviour: Antarctica, SRG, Pyrenees, bipolar south and
  REMBO `equil` 10 yr; ESM `equil` 1 yr; Greenland and bipolar north `none`;
  Laurentide and North `recon` (ICE-6G_C, regions 1.1/1.11/1.12). A transient
  Laurentide run (was "grow from zero ice") now sets `init_method = "recon_ref"`;
  a transient North run, `none`.
- The domain type `ice_domain` is renamed `kryos_domain`, in line with the
  Kryos naming of the cryosphere-component framework.
- `libs/yelmox_domain.f90` is split, by concept, into `kryos` (domain type,
  configuration, init, `remap`), `kryos_regions` (region-specific masks and
  physics), `kryos_coupling` (`step_*`, `couple_*_to_yelmo`), `kryos_startup`
  (cold start, restart bundles), `kryos_output` and `kryos_forcing` (`tsforcing`).
  Code is moved unchanged; drivers import each name explicitly.
- Renamed coupling primitives: `refresh_htopo` -> `couple_yelmo_to_htopo`,
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
  as one optional argument instead of `dTa`/`dTo`/`dSo`; the backend applies
  its anomalies only when it is active.
- `update_climate` is folded into `step_climate(dom, ts, tsf, init)`; `init=.true.`
  (the cold start) runs the update regardless of the `dt_clim` cadence.
- `domain_ctl` grid names: `grid_name` -> `grid_hub` (the hi-res hub),
  `grid_yelmo` -> `grid_ice` (Yelmo).
- Cold starts made consistent across drivers. `yelmox_esm` and `yelmox_rembo`
  now set up isostasy through the shared `domain_init_isostasy` (conservative
  ice-load coarsening + isostasy reference check). The optimisation's
  cold-start `cb_ref` (`domain_opt_init_cb_ref`) is set before
  `yelmo_init_state` in every driver (was after it in `yelmox`/`yelmox_bipolar`),
  so results of `opt` cold starts change slightly.
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
- Cold start: the surface is updated once more on the final initial geometry (after the
  `[sim]` ice init), so chion's active columns match the initial ice; before, new margin
  ice had smb = 0 until the next surface update.
- A forced run (`[tsforcing] active`) can restart from a bundle without tsgen state (an
  unforced spin-up): the forcing series starts at the restart time instead of stopping
  in ncio.
- `scripts/ant-paleo/run_lgp.sh` passes the spin-up restart as an absolute path
  (the relative one did not resolve from the run dir) and stops if it is missing.
- `[domain] regions_var`/`basins_var`/`sectors_var` left blank are empty (were
  undefined).
- Restarts continue the run exactly (bit-identical to the continuous run; needs
  FastIsostasy dev with restart-dt): the isostasy restart carries the ODE
  solver's time step and state (`ode_dt`, `ode_x`), and the sea level is restored
  as saved (`bsl_startup` no longer calls `bsl_update`; `A_ocean` restored).
  Before, the solver restarted from `dt_init`, giving tolerance-level differences
  that grew through the ice sheet.
- `esm_forcing`: in the historical period the direct-SMB anomaly is the ESM's SMB
  minus its own reference-period mean (`smb_esm_ref`), as for temperature,
  precipitation and the projection period. It subtracted the observed reference
  (`smb_ref`). Only transients with `use_esm`, `use_hist` and `use_smb` are
  affected; no current configuration runs that way.
- `yelmox_esm`: with `marine_shelf.extrap_shlf = True`, the reference ocean is
  extrapolated into the ice shelves on its own depth axis (`to_ref`/`so_ref`).
  It used the axis of the variability fields, which are loaded only with
  `use_var = True` (e.g. `scripts/1pctCO2/opt_ant.sh` sets `extrap_shlf` without it).
- `input/esm/esm_ant_ismip7.nml`: the SMB reference (`gcm_smb_ref`, read with
  `esm.use_smb = True`) is the RACMO2.3 monthly climatology `{grid}_RACMO23-VW23.nc`;
  the ERA5 1979-2022 file it pointed to has no `smb`.
- `opt.cf_init <= 0` starts the optimization from the till friction of the bed
  again (`cb_ref = cb_tgt`, from the `ytill` parameters), in every driver
  (`domain_opt_init_cb_ref`). Only `yelmox_rembo` still did; `yelmox`,
  `yelmox_bipolar` and `yelmox_esm` set `cb_ref = cf_init`, a negative friction.
  `cf_init > 0` is unchanged (uniform `cb_ref`).
- A restart from a bundle initializes Yelmo's passive-tracer backends (elsa,
  tracer): `domain_startup` loads Yelmo with `yelmo_restart_init` instead of
  `yelmo_restart_read`. Before, a restart with `ytrc.use_elsa` or `use_tracer`
  crashed in the first step. Requires a yelmo with `yelmo_restart_init`.
- `check_isostasy_reference` compares the two reference bedrocks on the isostasy
  grid, where Yelmo's `z_bed_ref` is remapped exactly as the isostasy reference was
  built: the same bedrock agrees to round-off on any isostasy grid (max |diff| <=
  1 m). Before, it remapped back to the Yelmo grid and allowed a 10 m mean
  difference, which a coarse isostasy grid over rough terrain exceeds (SRG on
  16 km: -14 m).
- SRG (Patagonia) runs again: `maps/grid_SRG-250M.txt` describes its grid (UTM
  zone 18S, 250 m), which the multigrid setup needs. Isostasy is off by
  default; when on, it runs on the new `SRG-16KM` grid (5x3 cells over the
  domain; at 250 m the padded FFT domain did not fit in memory). Requires
  fesm-utils with the transverse Mercator projection, and the
  `ice_data/SRG/SRG-250M` files with ascending `yc` (flipped on 2026-10-01).
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
- The original single-grid scripts and par files kept for the ports
  (`scripts/ant-paleo/legacy/`, `scripts/ismip7/legacy/`).

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
