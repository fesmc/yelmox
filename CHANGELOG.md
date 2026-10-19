# Changelog

All notable changes to YelmoX are recorded here. Each version corresponds to an
annotated git tag. Dates are release (tag) dates.

## [Unreleased]

### Added
- Prescribed ice `[ice_bnd]` (`libs/ice_bnd.f90`, optional): Yelmo's `H_ice_ref` in its
  fixed and relaxation cells (e.g. the yelmo `yelmo_masks.border_width` strip) from an ice
  history, `method = "file"` (at the forcing time, clamped to the file's range) or
  `"scaled"` (an index times one slice); `step_ice_bnd` in the drivers before
  `step_icesheet` and at the cold start (after the isostasy reference). `clim_db.nml`:
  `paleomist_H_ice` (PaleoMIST, 80-0 ka), `glac1d_H_ice` (GLAC-1D, 26-0 ka). Needs yelmo
  with `border_width`. `H_ice_ref` is in `yelmo.nc`. `input/yelmo_defaults.nml` synced.
- Climate modifier `[taper]` (optional), between `ocn_obm` and `bias`: the atmosphere
  towards the reference with a static weight field, `now = ref + w*(now - ref)` (the
  reference at the current surface); refused for rembo.
- `yelmox/yelmox_GreenlandPaleo.nml`: base config on GreenlandPaleo (GRL-PAL-16KM): ERA5 +
  ORAS5 present day with the MAR - ERA5 bias over Greenland (`[bias]`), "infinite" borders
  (`[yelmo] experiment = "infinite"`, `mask_border = "dynamic"`: the border ice follows its
  inner neighbours), sea level from file, uniform mantle. Prescribed border ice (a
  `border_width` strip with `[ice_bnd]`, PaleoMIST) is kept as an option, off.
- Greenland paleo, last glacial period (docs page `experiment-grl-paleo.md`):
  `yelmox_GreenlandPaleo_spinup.nml` (15 kyr, `cb_ref`/`tf_corr` optimization),
  `yelmox_GreenlandPaleo_lgp.nml` (-130 kyr to 0, PMIP3 piControl/LGM on the glacial index
  `alpha_combined_125kyr_interp`, `ocn_fraction`, transient insolation, elsa),
  `scripts/grl-paleo/run_spinup.sh`, `run_lgp.sh`.
- MAR v3.14 - ERA5 1981-2010 bias on GRL-PAL for `[bias]`: `scripts/grl-paleo/make_bias_mar_era5.jl`
  (difference at sea level with the `[snap]` lapse rates and `f_p`; outside the MAR coverage
  the 50 km-smoothed bias tapered to zero over 100 km), written to
  `ice_data/GreenlandPaleo/GRL-PAL-16KM/`; `clim_db.nml` groups `bias_mar_era5_1981_2010_dtas`,
  `_fpr`.
- `input/clim_db.nml`: shared climate database of the snapshot backend, one varslice group per
  product field. Present-day references regional (`{domain}/{grid_name}`): MAR v3.14
  1961-1990, RACMO2.3p2 1981-2010, ERA5 1981-2010, ORAS5 1981-2010, ISMIP7 ocean.
  Paleoclimate from the global FesmData products, remapped online (`bilin`): PMIP3 and
  PMIP4 ensemble means, CLIMBER-3α, CLIMBER-X (with the anomalies to piControl), atmosphere
  and ocean. Path keywords `{grid_global}` (`[snap] grid_global`, default `GLOBAL-1DEG`)
  and `{exp}` (the snapshot's `exp`).
- `input/yelmox_defaults.nml`: defaults of the forcing groups, `[indices]`, `[idx_<name>]`
  (template `[idx]`, passed to `tsgen_init`), `[snap]`, `[snap_field_<name>]`,
  `[snap_<snapshot>]`, `[snap_axis_<name>]` (templates `[snap_field]`, `[snap_snapshot]`,
  `[snap_axis]`), `[anom]`, `[ocn_fraction]`, `[bias]` (domain suffixes on the template);
  read with `nml_read(..., defaults_file, defaults_group)`, unknown keys stop the run
  (`nml_validate`). Needs fesm-utils `tsgen_init(..., defaults_file, defaults_group)` and
  `nml_has_group` (fesm-utils `forcing`, to be merged into dev).
- `scripts/check_par_nml.py` also checks the forcing groups against
  `input/yelmox_defaults.nml` (named groups against their templates) and the snapshot
  bindings against the `var_defs` files; `--dead` reports unread snapshot groups.
- `scripts/climate-compare/`: climate-only comparison of two par files (`run.sh`: ice sheet
  off, `tstep_method = "rel"`, climate every 1 kyr over -130 ka to 0 or PD; `compare.jl`:
  domain means per time and maps of `t2m_ann`, `pr_ann`, `to_ann`), the regression check
  for climate backend and modifier changes.
- Static climate bias correction `[bias]` (optional): additive temperature and
  multiplicative precipitation fields (monthly or annual, varslice groups, online remap
  allowed) with weights `f_tas`, `f_pr`, applied to the monthly atmosphere of snapclim,
  snapshots and esm before downscaling; written to the climate 2D file.
- The snapshots backend reads its snapshots on the climate grid (`snapshots_init` takes the `grid_class`),
  so `var_defs` groups can be remapped online (`remap`, `grid_src`).

### Changed
- Snapshots: one monthly lapse rate both ways. The inflation of the blended sea-level
  fields to the model surface (`snapshots_transform`, also the surface hosts' monthly lapse
  rates) used the seasonal cosine at phase `m*30-15`, the sea-level reduction of a snapshot
  at `m*30-30` (as the season synthesis, `cos(2π(m-1)/12)`: summer lapse rate in July, or
  January in the south). Both now use `lapse_monthly` on the calendar of the synthesis, so a
  snapshot at its own elevation is recovered exactly (the `taper` reference was off by
  `z_srf*(lapse(1)-lapse(2))*(cos(φ-15°)-cos(φ-30°))`). Monthly `tas` at the surface changes by
  up to ~0.4 K per km of elevation in spring/autumn (sign alternating, zero in the annual
  mean to first order); summer `tas` and `smb` change slightly.
- Index window modifier (fesm-utils tsgen `win_t0`, `win_t1`, `win_taper`, `win_fac`, `win_ref`,
  `win_above` in `[idx_<name>]`; replaces the snapclim `f_hol` rule) and the RNG state in
  the index restart files (noisy indices continue identically after a restart). Needs fesm-utils
  with these changes.
- Insolation from the [insol](https://github.com/fesmc/insol) library, now a dependency
  (`.configme/manifest.toml` `insol:main`, `make insol-static`, `libinsol.a`); the vendored
  copy `libs/insol/` (and its `interp1D`, smbpal's `interp_time` now uses fesm-utils') is
  removed. The orbit is `[orbit<sfx>]` (`libs/orbit_forcing.f90`, defaults in
  `input/yelmox_defaults.nml`), shared by chion's host and smbpal: `method = "clock"` (at the
  forcing clock), `"time"` (at a fixed `time`) or `"fixed"` (`ecc`, `obl`, `per`); it
  replaces `const_insol`, `const_kabp` and `insol_fldr` of `[surface_chion]` and `[smbpal]`
  (`const_insol = True, const_kabp = X` -> `method = "time", time = -X*1e3`;
  `const_insol = False` -> `"clock"`). smbpal computes its insolation on the grid and only
  when the orbit changes. The insolation changes by up to ~2.5 W/m², mostly at low
  latitudes: insol integrates the daily mean exactly, the vendored copy summed 24 hourly
  samples. Existing builds: `make clean` once (stale `insolation.mod`, `interp1D.mod`).
- `yelmox/run_mg_forward.sh`: the anomaly from `[anom]` on the constant index `[idx_dT]`
  (`yelmox_Antarctica.nml`) instead of the removed `snap.dTa_const`/`dTo_const`.
- snapclim removed (`libs/snapclim.f90`, `[comps] climate = "snapclim"`, its tests,
  `input/greenland_snp.nml`, `yelmox_Greenland_snapclim.nml`); the snapshot backend snapesm
  renamed `snapshots` (`[comps] climate = "snapshots"`, `libs/snapshots.f90`,
  `snapshots_*`; groups keep `[snap*]`, output `snap.nc`). `[snap_field_<name>] kind`
  removed (the field name decides). `yelmox_Greenland_snapesm.nml` removed (its climate is
  that of `yelmox_Greenland.nml`), `yelmox_Greenland_snapesm_ramp.nml` renamed
  `yelmox_Greenland_ramp.nml`. Docs page `climate-snapshots` (was `climate-snap`, redirected).
- The regional configs run on snapshots (`[comps] climate = "snapshots"`, snapclim groups
  removed): `yelmox_pd_Greenland`, `yelmox_pd_Antarctica`, `yelmox_Antarctica` and
  `yelmox_Antarctica_paleo_spinup` (`pd` snapshot only; spin-up with `[ocn_fraction]`
  f_to 0.40), `yelmox_Antarctica_paleo_lgp`, `yelmox_North`, `yelmox_LIS` (`pd`,
  `picontrol`, `lgm` on an axis `a`: LGP index `alpha_combined_125kyr_interp` as is,
  North/LIS `const` 1; `[ocn_fraction]` f_to 0.25), `yelmox_Greenland` (as
  the former `yelmox_Greenland_snapesm`), `yelmox_Pyrenees` (`pd` + `[anom]`, index
  `alpha_orbital_t2018` x f_ta 10 K), `yelmox_SRG` (`pd`, the legacy MAR v3.14 file as
  groups of the par file). Physical changes: Antarctic reference from FesmData RACMO2.3p2
  1981-2010 t2m/pr (was the RACMO2.3/ERA-Interim hybrid: -1.0 K, +2.5% precipitation
  over the ice) and the ISMIP7 ocean climatology (was ISMIP6 J20: +0.1 K); PMIP3 from
  the global product (GRL +0.3 K over the ice); the ocean anomaly from the sea-level air
  temperature anomaly without lag (North/LIS LGM ocean 0.8/1.8 K warmer: snapclim's rule
  included the lapse-rate cooling of the ice sheets); dropped: `f_hol` (GRL), the 80 km
  smoothing of the reference before 1 ka (`f_corr`: LGP, North, LIS; local +-4 K, means
  unchanged), the unused snapshot oceans, the SRG and Pyrenees oceans (`with_shelf =
  False`).
- yelmox_bipolar: the OBM box air temperatures from one global anomaly index,
  `[obm_forcing]` (`index`, `gain`): `theta_i = theta_i_init + ampl_i × gain × index` for
  every box (the tropical box set once), vapour fluxes likewise; `climate_air_anom` and
  `libs/obm/atm2ocean.f90` removed. `obm2ism` is the climate modifier `[ocn_obm_<sfx>]`
  (`absolute`, as before, or `anomaly`: reference + box anomaly), applied by the driver;
  `[ctrl] obm2ism` removed. The OBM restart keeps the `*_init` values of the par file and
  now holds the box air temperatures. The Bipolar par file on snapshots (north as GRL, south
  present day), indices `at` and `gmt`.
- `[snap] var_defs` is a list of files, searched in order (first that has the group; an
  error names the group and the files), then the par file itself (implicitly last, so a
  par file never names itself); also `[bias] var_defs` (default `""`: the par file only). The forcing groups need only the keys that
  differ from `input/yelmox_defaults.nml`: a field without axes needs no group, `[anom]`/`[ocn_fraction]`/
  `[bias]` absent = inactive. GRL configs take PMIP3 piControl/lgm from the global
  product (`exp`, remapped; regional MAR/ORAS5 unchanged); `input/greenland_clim.nml`
  removed.
- snapshots, phase 3: fields optional per snapshot (key absent or `""`); the backend builds
  what the reference binds (atmosphere `tas`+`pr`+`zs` and/or ocean `to`, `so` optional),
  bindings checked at init (axis snapshots bind the fields they blend, `zs` with `tas`/`pr`,
  group counts); ocean annual (one slice, e.g. CLIMBER ocean products) or a monthly
  climatology, salinity 0 without `so`; axes with `lo = ""` (`hi` is the anomaly or ratio,
  without `zs` taken at sea level); `tas_season = "jja" | "djf"` for two-field `tas`
  (default the hemisphere's summer; PMIP3 `clim_*_tas_djf` groups added); `ref_name` must
  name a snapshot; `[snap] f_stdev` removed; `snapshots_end` and `climate_end` (drivers).
  rembo takes its ocean from snapshots, an ocean-only `pd` snapshot (ORAS5; the shipped
  config had no ocean, `ocn_path = "none"`), plus `dT_ocn`; `climate = "snapshots"` needs an
  atmosphere.
- Climate modifiers, phase 2: `[ocn_fraction<sfx>]` (`active`, `f_to`, `region`) after
  `[anom]`: `to_ann = to_ann_ref + f_to × mean(tsl_ann − tsl_ann_ref)` over `region`
  (`"all"` or a `[regions]` selection, resolved at init), without the lag of one climate
  update and on the sea-level anomaly (no lapse-rate cooling of the ice sheet); refused
  for esm and rembo and with `[anom] f_to /= 0`. It replaces the snapshot backend's `blend = "fraction"`
  (now an error; `to` with `blend = "const"`), `[snap] f_to` and the lagged `ta_ann`,
  which changes the results of every `fraction` config (`yelmox_Greenland*.nml`).
  The bias correction is the last modifier, group `[bias]` (was `[climate_bias]`, same
  keys); `anom` with `f_to`/`f_so` /= 0 on a climate without ocean stops at the first
  update. Climate 2D file: `t2m_sum`, the summer mean of the applied `tas`, for every
  backend with monthly fields (esm: was the reference plus `dts`; rembo keeps REMBO's);
  snapshots no longer writes `tsl` and `ta_sum`; `ocn_fraction_dto`.
- Climate modifiers (`libs/climate_mods.f90`): optional groups applied to the backend's
  output in a fixed order, contracts checked at init. First one `[anom<sfx>]` (`active`,
  `index`, `f_ta`, `f_to`, `f_so`; replaces `[tsforcing]`): homogeneous dTa/dTo/dSo after
  the backend (snapshots: no longer before the elevation transform, equal up to rounding;
  its fraction ocean no longer sees dTa), before `[bias]`; refused with the esm
  smb. rembo reads `[anom]` itself. The recon cold start (`init_method = "recon"`) is now
  forced too. Removed `[tsforcing]`, the dTa/dTo/dSo arguments of snapshots and its
  `[snap] dTa_const/dTo_const/dSo_const`. snap.nc `tas`, `pr`, `tsl_ann` are the applied
  climate (after the modifiers).
- Indices registry (`libs/kryos_forcing.f90`): named scalar series `[indices] names`, each a
  tsgen group `[idx_<name>]` (optional `units`, `restart_every_df`), owned by the driver,
  advanced once per step on the forcing clock `time_rel` with the ice volume as response,
  shared by the domains of `yelmox_bipolar`; kill switch and Δf restarts per index; restart
  files `idx_<name>_restart.nc`, 1D output `idx_<name>`. snapshots axes `[snap_axis_<name>]`
  (`lo`, `hi`, `index`), listed in `[snap] axes`; fields list `axes` (was `index`); climate
  output `axis_<name>` (was `idx_<name>`). `[tsforcing]` takes its value from an index
  (`index`). Removed `[tsgen]`, `[snap_idx_*]`, `[tsforcing] restart_every_df`,
  `tsgen_restart.nc` and `snapesm_idx_<name>_restart.nc`. The ramp config runs with
  `tstep_method = "rel"` (the ramp now advances on `time_rel`).
- The climate 2D file has the annual ocean profiles `to_ann`, `so_ann` (over `depth`) for
  every backend that has them (snapclim, snapshots, rembo); written by `climate_write_2D`
  (were snapshots only).
- Climate on a coarser grid than the surface: `tas`/`pr` are moved from the climate grid's
  surface elevation to the surface grid's with the backend's monthly lapse rate and
  precipitation sensitivity (`climate_out`: `now%z_srf`, `lapse(12)`, `now%beta_p`;
  `kryos_coupling` `remap_atm_to_surface`). Unchanged when `grid_clim = grid_surface`.
- The hub geometry reaches `grid_clim` and `grid_surface` conservatively when the hub is
  finer (`remap_method_smooth`; was bilinear), which changes runs with a finer hub.
- snapshots indices are coordinates: each `[snap_idx_<name>]` spans an axis `snap_lo` ->
  `snap_hi` and its tsgen value is the position on it (any tsgen method, e.g. `const` for a
  steady state); a field may take several indices (`index = "at" "bt"`), whose axes are
  summed (`linear`) or multiplied (`ratio`) on the reference, as snapclim `snap_2ind`
  (Banderas et al., 2018). Normalising a proxy series is done in tsgen (`norm_t0`,
  `norm_t1`; fesm-utils). Removed `manifold`, `f_hol` and the snapshot keys `time`,
  `idx_coord`. Needs fesm-utils with tsgen normalisation.
- FesmData v2 products from the link `fesmdata` (`$FESMDATA`, e.g.
  `/work/ba1442/fesmdata`; was `ice_data/v2`), copied into the run folders by
  runme (`.runme/info.json` links). SRG topography `TOPO-GLO90-GEBCO2026`.
- Isochrones from `Stratigraphy-MacGregor2025` (`pd_age_names = "age" "depth_iso"`; was
  `STRAT-M15`). `input/yelmo_defaults.nml` synced with yelmo.
- Pyrenees on FesmData: ice sheet and surface on the hub `PYR-250M` (`TOPO-GLO90-GEBCO2026`),
  climate (`Atmosphere-ERA5-1981-2010`, was WorldClim) and isostasy on `PYR-2KM`, GHF
  `GHF-Lucazeau2019` from `PYR-2KM`; no ocean, no pd T_srf/velocity (the PYR-1KM `ice_data`
  folders are empty). Climate and surface every 50 yr (`dt_clim`), chion `time_equil = 10`
  (chion costs ~15 s per year on PYR-250M).
- GHF and sediment obs are remapped (bilinear) onto the Yelmo grid when the file is on another
  grid (FesmData global attribute `grid_name`; fesm-utils `fesmdata_read_field`), so a par file
  can name a coarser data grid in `obs_path`. `geothermal_init`/`sediments_init` take the Yelmo grid.
- Greenland present-day climate MAR v3.14 1961-1990 (FesmData `SMB-MARv3.14-1961-1990`, was
  MARv3.11-ERA): snapclim/snapshots reference climatology (`t2m`, `sf`, `rf`, `z_srf`), the esm
  reference blocks (`esm_grl_*.nml`) and the yelmo pd smb/T_srf (`-annual-mean`, `T_srf` skin
  temperature). snapclim reads precipitation and temperature units from the files.
- FesmData v2 velocities (`Velocity-Joughin2018`, `Velocity-Rignot2011`), GHF
  (`GHF-Shapiro2004`, `-FoxMaule2005`, `-Loesing2021`, `-HazzardRichards2024`; error `ghf_sd`)
  and North/LIS reference climate `Atmosphere-ERA5-1981-2010` (was ERA-Interim). rembo1 on
  `TOPO-BedMachine-v6` and MAR v3.14. Unused esm `&basins` (BASINS-nasa) removed.
- Paleo climate snapshots from FesmData `Paleoclimate-PMIP3-{piControl,lgm}-{mean,stdev}` (was
  `PMIP3_sig50km`, smoothed with 50 km) and ocean from `Ocean-ORAS5-1981-2010` (was ORAS4): snapclim,
  snapshots (`greenland_clim.nml`) and the esm ocean reference (`esm_grl_*.nml`, monthly climatology,
  `fill = "nearest"`; read from `GRL-16KM` and remapped bilinearly onto the run grid, since ORAS5
  is cropped on GRL-16KM/32KM only). snapclim reads ocean fields with missing values and fills them
  (fesm-utils `fill_nearest_levels`: each level from the nearest valid cell, else from the level
  above); the ocean mask comes from the data where `ocn_names(2) = "none"` (PMIP3 `sst_ann`).
  A 2D ocean without salinity gets 0 (salinity was read into an unallocated array).
- snapclim/snapshots snapshot keys `clim_sigma` and `ocn_sigma` [km] (optional, default 0): Gaussian
  smoothing (fesm-utils `filter_gaussian`) of the atmosphere / ocean of the snapshot at loading.
  The PMIP3 snapshots use 50 km (as the v1 `PMIP3_sig50km` files); `snapclim_init`/`snapshots_init`
  take the grid spacing.
- chion replaces smbpal in the PDD configurations (`yelmox_Antarctica`, `yelmox_pd_Antarctica`,
  `yelmox_North`, `yelmox_Pyrenees`, `yelmox_esm_Antarctica_nudge`, bipolar south):
  `surface_method = "chion"`, `[chion] model = "pdd"`, `[pdd_chion]` from smbpal's PDD
  (`pdd_method = "cg05"`, Calov and Greve, 2005). The chion host supports PDD (no
  insolation), and takes the surface temperature of ITM and PDD from chion's
  `chion_tsrf_annual` (smbpal's firn warming by the annual refreezing, on annual means).
  Needs chion with PDD melt/refreezing, `cg05` and `chion_tsrf_annual`;
  `input/chion_defaults.nml` synced.
- Surface output: one file per surface model, `<surface_method>.nc` (`chion.nc`,
  `smbpal.nc`, `smb_simple.nc`, new); with `surface_method = "climate"` none, the climate's
  file has the smb applied (`smb_ann`, `tsrf_ann`, for esm and rembo). The climate's smb is
  landed on Yelmo from the climate (no longer through smbpal), and smbpal is initialized
  and restarted (`smbpal_restart.nc`) only when it is the surface model.
- ISMIP7 (`yelmox_esm_Antarctica_ismip7.nml`, `yelmox_esm_Greenland.nml`,
  `yelmox_esm_Greenland_k24.nml`): `[comps] surface_method = "climate"`, the ESM smb of
  the protocol (was chion since the chion default, smbpal before), so that
  `run_ismip7_*.sh` and the `opt_*.sh` spin-ups (which set it themselves before) impose it.
- ISMIP7 Antarctica: the reference climate of `esm` (`t2m`, `pr`, `smb`, `z_srf`) from
  FesmData v2 `SMB-RACMO2.3p2-1985-2014`, the RACMO2.3p2 monthly climatology of the ISMIP7
  reference period, complete over the grid (was the RACMO2.3p2 ERA5 1979-2022 series and
  the VW23 climatology).
- Antarctica: the observed basal melt (`[marine_shelf] use_obs`) from FesmData v2
  `BasalMelt-Rignot2013`, `bmelt_eq_reese` (was `BMELT-R13`, `bm_eq_reese`), missing
  values filled from the nearest valid value; Greenland, North, LIS and SRG keep `BMELT-R13`.
- Sediment thickness (`[sed] obs_path`) from FesmData v2 `Sediment-CRUST1` (CRUST1.0,
  Laske et al., 2013, successor of the Laske and Masters (1997) `SED-L97`); SRG keeps
  `SED-L97` (no v2 grid).
- ISMIP7 Antarctica: present-day smb and surface temperature from FesmData v2 RACMO2.3p2
  1981-2010 (`ice_data/v2/`); `T_srf` is now the skin temperature (was 2 m). The smb is in
  kg m-2 yr-1: needs Yelmo with the pd smb units read from the file.
- **FesmData v2 regions** (needs FesmData v2 data, `ice_data/v2/`): the topography is in
  `[topo]` and the regions, zones and basins in `[regions]` (fesm-utils `topodata`,
  `regions`); htopo moved to fesm-utils. Each component makes its own masks from the
  regions on its grid with selection expressions: Yelmo in `[yelmo_masks]` (where ice is
  allowed, relaxation, error metrics, basins, named regions), the marine shelf in
  `[marine_shelf]` (`mask_ocean`, `mask_deep_ocean`, `mask_pico_deep`, `basins`, and the
  corrections `corr_method`/`corr_names`/`corr_values`/`corr_<name>` replacing `tf-grl`,
  `tf-ant`, `basin_number`, `basin_*_corr`), the climate in `[esm] basins`, the domain in
  `[sim] recon_regions` and `[negis] region_*`. Removed from `[domain]`: `topo_*`,
  `regions_*`, `basins_*`, `sectors_*`, `ice_codes*`, `region_*`, `relax_*`. The
  `tf-ant` corrections keep the Zwally ids they read before (`Zwally2012:1,12,14,15`),
  although their labels follow the Reese numbering.
- Existing builds need `make clean` once after updating: the `htopo.mod` left in
  `libyelmox/include` by the old yelmox htopo hides the fesm-utils module, and `kryos.f90`
  fails to compile (`'ref' is not a member of the 'htopo_class' structure`).
- Removed the unused NE precipitation factor of snapclim (`f_p_ne`) and its basin masks.
- Par files: `ytopo.grad_lim_zb = 0.5` (Yelmo default; was 0.05, except SRG 0.5 and
  Pyrenees 1.0). With 0.05 Yelmo smoothed the bed and ice thickness wherever the bed slope
  exceeds 5% (initial state and present-day target): ANT-16KM grounded thickness 131 m RMSE
  from the observations at the start of the ISMIP7 spin-up; with 0.5, 16-24 m (56-59 m after
  15 kyr, was 141 m). Ice-sheet results change.
- esm: the ocean fields are extrapolated into the ice-shelf interiors on copies, not in
  the forcing slices (fesm-utils varslice now reads data without time once). Fixes the
  repeated extrapolation of slices that were not re-read (ESM reference ocean of the
  transient runs, variability) with `extrap_shlf = True`; results otherwise unchanged.
- Model time in double precision (needs fesm-utils, Yelmo, FastIsostasy and chion with
  `real(dp)` time). The drivers' step `dtt`, `[comps] dt_clim`/`dt_isos`, `[output]
  dt_cmip` and the time arguments of the kryos modules and of the climate, surface and
  shelf components (snapclim, snapshots, esm, rembo, smbpal, chion surface, marine_shelf,
  OBM coupling) are `real(dp)`, as are the snapclim series time axes and snapshot times
  and the esm reference periods. Forcing values stay `wp`. In single precision the
  spacing of the time is 0.008 yr at 130 kyr.
- `input/yelmo_defaults.nml` synced with Yelmo dev (`[ytrc]` horizontal advection keys,
  needed to run; `bkt_floating_mode` default 2, the par files keep 0).

### Added
- `yelmox_bipolar`: the ocean box model alone (`active_obm` without `active_north`/
  `active_south`); `atm2obm`, `ism2obm` and `obm2ism` then stop the run (they need an
  ice-sheet domain).
- `[comps] dt_isos` [yr]: isostasy update interval (10 in all par files; was every step).
  FastIsostasy integrates the interval in one call; the bedrock is held in between. Its
  yearly `tsit54` step was ~37% of the ANT-16KM ISMIP7 spin-up wall time (11% at 8KM).
  Results change slightly.
- ESM spin-up scripts (`scripts/ismip7`, `1pctCO2`, `tipmip`) pass `comps.dt_clim=10`
  (the forcing is a fixed climatology; ~12% of the ANT-16KM spin-up wall time at 1 yr);
  transient runs keep `dt_clim = 1`.
- Timing of the time loop: the wall time of each part (output, forcing, relax, optimize,
  isostasy, coupling, ice sheet, climate, surface, shelf; bipolar also the ocean box model)
  is printed at the end of a run and written to `timing.txt` at the `tm_1D` cadence
  (accumulated seconds). Needs fesm-utils with `timer_init`/`timer_write_totals`.
- chion BESSI in yelmox (`[chion] model = "bessi"`): the host supplies its surface shortwave
  (TOA insolation times `[surface_chion] trans_sw`), constant `wind_speed` and `rel_hum`, the
  air pressure from the surface elevation (barometric, annual-mean air temperature) and the
  solar longitude; longwave from chion's own parameterization. ITM results unchanged.
- 1pctCO2 (Greenland, Antarctica) and TIPMIP (Greenland) finished: ice sheet and isostasy
  on, ISMIP7-like optimization spin-up; TIPMIP stabilisations branch off the ramp's
  restarts (years 109, 232). 1pctCO2 reads the global CMIP fields of `[esm] esm_name`
  (default MPI-ESM1-2-LR) with online remapping (fesm-utils `varslice` `remap`/`grid_src`).
- The geothermal heat flux (and its error field) is converted to mW/m2 from the `units`
  attribute of its variable (W/m2 or mW/m2; other units stop the run), so the same par
  file reads either version of a dataset; missing (negative) GHF values are filled from
  the nearest valid value.
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
- snapshots writes its state to `snap.nc` (`[output] write_clim`): the driving
  indices `idx_<name>`, `z_srf`, monthly `tas`/`tsl`/`pr`, `ta_sum`, `tsl_ann`
  and the ocean profiles `to_ann`/`so_ann` (before only `t2m_ann`, `pr_ann`).
- snapshots restarts: restart bundles hold the state of each driving index
  (`snapesm_idx_<name>_restart.nc`, with the snapshots configuration as
  attributes), restored on restart; bundles without them keep the cold-start
  indices.
- `scripts/ismip7/`: the ISMIP7 optimization spin-ups (`opt_ant.sh`, `opt_grl.sh`,
  and L. Gutierrez Gonzalez's `opt_grl_ismip.sh`) for the current `yelmox_esm`.

### Changed
- chion BESSI follows chion `sync-dev-nils` (Chion.jl 9ec6cc7; `input/chion_defaults.nml`
  re-synced): semix turbulence, cloud-proxy longwave, 5-layer ice substrate, thin-snow albedo
  blend, 8 diurnal substeps. yelmox now fills `H_ice` for BESSI (`H_ice = 0` is land: no ice
  albedo, substrate or ablation) and passes its TOA insolation to the cloud proxy. Greenland
  par files: `[chion_const] alpha_ice = 0.55` (BESSI; chion default 0.40), calibrated to MAR
  v3.11 1961-1990 like ITM's `alb_ice`: runoff 648 Gt/yr (MAR 645; 0.50 gave 705), ablation
  area 0.43 (MAR 0.41) 1e6 km2, ablation-zone smb bias +15 mm/yr. ITM results unchanged
  (bit-identical). BESSI GRL-8KM (16 threads): cold-start setup 664 s (was 212), 20-yr loop 2.1x.
- Greenland par files: chion's ITM parameters in their own group `[itm_chion]` (`[chion]
  nml_itm`), `[itm]` with `alb_ice = 0.31` (was 0.4, shared with smbpal), calibrated to MAR
  v3.11 1961-1990 at GRL-8KM with fixed present-day geometry: runoff 635 Gt/yr (MAR 645;
  0.4 gave 543), ablation-zone smb bias -2 mm/yr (0.4: +215). smbpal keeps `[itm]`.
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
- Run scripts (`scripts/ismip7`, `ant-paleo`, `1pctCO2`, `tipmip`, `run_lis_smbsimple.sh`): the
  queue options are taken from the environment if set (`runopts`, `HPCOPT_SPINUP`/`HPCOPT_SCEN`,
  `RUNOPT`), so they submit on other machines than Levante (e.g. albedo `-q 48h`).
- `scripts/check_par_nml.py` follows yelmo's defaults file: it reports keys of yelmo groups
  missing from `input/yelmo_defaults.nml` (what `nml_validate` stops on) and duplicate
  groups/keys, no longer omitted parameters; `--dead` also reports non-yelmo keys that no
  source reads. Checks `yelmox/*.nml` and `yelmox_bipolar/*.nml` by default.
- **ISMIP7 geothermal heat flux**: the ObsISMIP7 GHF (W/m2) was never converted and,
  with `obs_err_name` = the GHF itself, clamped to 0.1 mW/m2 everywhere. Fixed
  (units from the file, `obs_err_name = "none"`); ISMIP7 spin-ups must be redone.
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
- The hemisphere of a domain (seasons, lapse rates in snapclim, snapshots and the
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
  (was `write_snap`): `snap.nc` (snapclim, snapshots) or `esm.nc` (the esm fields:
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
  "snapshots"` replaces the `make CLIMATE=` switch (one `libs/yelmox_climate.f90`
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
- `maps/grid_SRG-250M.txt` is the FesmData v2 SRG-250M grid (208 x 120, as the product
  regenerated on 2026-10-16, which now shares the grid of the legacy `ice_data` SRG inputs);
  it held a stale 257 x 129 description. `yelmox_SRG.nml` runs every component on
  `SRG-250M`; `maps/grid_SRG-250M-v1.txt` removed.
- Exact restarts (bit-for-bit with the continuous run, `yelmox` and `yelmox_bipolar`): the
  restart path re-ran `step_surface`, stepping the chion/smbpal snowpack one extra year
  (smb, refreezing, H_snow off by up to ~100 %). The chion and smbpal restart files now
  carry the annual fields of the last cycle (`ann_*`), and `domain_startup` continues the
  restored state, re-deriving only the climate (and a stateless `smb_simple` surface); the
  drivers no longer re-step surface and shelf (`obm_ocean_init` steps the shelf on a cold
  start only). Bundles without `ann_*` still restart, inexactly, with a warning. Needs
  FastIsostasy dev ca616b3 (each domain keeps its own restored bsl, which the
  continuous multi-domain run sees). Remaining differences at the restart record only: the
  Yelmo diagnostics `visc_bar`, `mb_clip_tot` (not in the Yelmo restart).
- `const_kabp` (surface_chion, smbpal) is kyr BP: the constant insolation used the orbit at
  +const_kabp kyr (future) instead of -const_kabp kyr. No effect for `const_kabp = 0` or
  model-time insolation.
- `tests/test_timestepping.f90`: current fesm-utils API (double-precision time, methods
  rel/cal/const); checks the end time, step count and calendar/relative times.
- `maps/grid_GRL-32KM.txt`: 53 columns, as the GRL-32KM data and the FesmData v2 grid (was 54;
  the ncio size check stopped `yelmox_bipolar` on the GRL-32KM climate). Delete cached
  `maps/scrip-*GRL-32KM*.nc` maps.
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
- `libs/ismip6.f90` (module `ismip6`, unused since the `yelmox_ismip6` and
  `yelmox_nahosmip` flavors were retired).
- `libs/hyster.f90` (module `hyster`, unused since tsgen replaced it) and the dead
  `yelmox_hysteresis_help` make rule (its source was never in the repository).
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
