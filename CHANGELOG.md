# Changelog

All notable changes to YelmoX are recorded here. Each version corresponds to an
annotated git tag. Dates are release (tag) dates.

## [Unreleased]

### Added
- Coefficients of kind `scaled` (`[coef_<name>] method = "scaled"`, `coef`, `factor`):
  `factor` x a tsgen coefficient of the registry, set right after it each step (terms on
  the atmosphere may use it), recomputed on a restart; `coefs_varies` and
  `coefs_check_used` in `libs/kryos_forcing.f90` (the coefficient kind is now
  `tsgen | scaled | derived`). `tests/test_coefs.f90` covers them.
- `yelmox_bipolar/params.sh`: the coupling setups of yelmox_bipolar as runme `-p` arrays
  (uncoupled, box model alone, box model -> both ice sheets, Greenland -> box model,
  Greenland <-> box model, Antarctica <-> box model, fully coupled) and the modifier
  `bip_ocean_anomaly` (the reference ocean plus the box anomaly). `run_tests.sh` runs
  every setup and the consistency checks (each must stop with its error), or the setups
  named (`setup+modifier`).
- Docs `yelmox-bipolar.md`: "Coupling setups", the exchanges (diagram), their switches,
  the setups table and the rules behind it.

### Changed
- Coefficient output: the driver writes the coefficients and the barystatic sea level to
  the run's own 1D file `forcing_ts.nc` (`[tm_1D]` cadence, both drivers; in
  yelmox_bipolar the run folder, not each domain's), under their names (`gmt`, not
  `coef_gmt`), with `kind` and the source (`method`, `coef` + `factor`, `from`) as
  attributes. `yelmo_ts.nc` no longer holds them. `time` and `bsl` are not coefficient
  names (checked at init).
- yelmox_bipolar: the regional amplification of `gmt` is set once, in the scaled library
  coefficients `dTa_north` (1.07), `dTa_tropics` (0.79) and `dTa_south` (0.86 x `gmt`;
  Zickfeld et al., 2004), read by the domains (`term_regional_north`/`_south` on them,
  their `f_ta` dropped) and by the box model (`[obm_forcing] coef_north`, `coef_tropics`,
  `coef_south`; `coef` the global anomaly of the vapour fluxes; `gain` removed);
  `[nautilus] thermal_ampl_*` removed. No change in answers.
- yelmox_bipolar: `[ctrl] ism2obm` removed; the freshwater flux of each hemisphere is
  switched by `couple_fwf_north`/`couple_fwf_south` alone (an error without the box
  model or without the domain of that hemisphere).
- Consistency checks at init: every coefficient of `[forcing] coefs` is read by the run
  (a term of an active domain, `[ice_bnd]`, the OBM forcing, or the source of a scaled
  coefficient that is; `coefs_check_used`, both drivers); in yelmox_bipolar, with the box
  model on, `atm2obm` when the atmosphere of a domain changes in time. yelmox_bipolar
  runs its checks before the startup of the domains (`setup_domain`, checks,
  `start_domain`).
- yelmox_bipolar: consistency checks of the ocean terms at init (`obm_coupling_init`): an
  error for a box coefficient of the other hemisphere (`obm_south` in the north), for a
  domain with its ocean both from a box and from the air anomaly (`obm_north` with
  `ocean`), and for an air-anomaly coefficient used by both domains (`dTo`: the
  coefficients are shared, the second domain would overwrite the first).
- yelmox_bipolar par file: the north without the term `lgm` and the coefficient `alpha`
  (zero at present day, so no change in answers): one coefficient, `gmt`, drives the box
  model and both domains (`[forcing] coefs = "gmt" "dTo"`, `[forcing_north] terms =
  "regional_north" "ocean"`); the forcing and coupling comments point to the setups.

## [v2.4] - 2026-10-23

### Added
- The forcing as terms (docs `yelmox.md`, `climate-snapshots.md`): the climate at the
  forcing time is a reference snapshot plus terms, each a coefficient (a scalar of the
  forcing time) times a pattern; the fixed numbers of the par file (the gains, beta_p, the
  lapse rates) are its parameters. `[forcing]`: `reference`, `coefs`, `terms`, one
  `var_defs` list for every named group (default `"input/clim_db.nml"
  "input/forcing_db.nml"`). `[coef_<name>]` (template `[coef]`): a tsgen coefficient
  (`method = "series"` a time series from a file, normalised with `norm_t0`/`norm_t1`;
  `const`, ramps, controllers; window `win_*`; Δf restarts `restart_every_df`), or `method
  = "derived"`, set by the model each update (`from = "air_anomaly"`: `factor` x the mean
  sea-level air anomaly over `region`; `"obm_north"`/`"obm_south"`: a box temperature of
  the bipolar ocean box model, its anomaly with `anomaly = True`); ocean fields only.
  `[term_<name>]` (template `[term]`): `coef` (`""` = 1), `pattern` (`snapshots`: `snap0`
  -> `snap1`, or `snap1` alone as a given anomaly or ratio; `uniform`: `f_ta`, `f_to`,
  `f_so`; `replace`: the ocean field = the coefficient), explicit `fields`; the rule fixed
  by field (tas, to, so add; pr a ratio, exp(beta_p f_ta f) in a uniform term). Contracts
  checked at init (coefficients and snapshots exist, fields bound, replace once per field,
  derived coefficients on the ocean, uniform atmosphere terms need an atmosphere that
  drives the surface or rembo, refused with esm and its smb; ocean terms an ocean). The
  coefficients are owned by the driver (`libs/kryos_forcing.f90`), advanced once per step
  on the forcing time `time_rel` with the ice volume as response, shared by the domains of
  `yelmox_bipolar`, evaluated once per climate update (`libs/climate_terms.f90`). Output
  `term_<name>` (2D), `coef_<name>` (1D, in `yelmo_ts.nc`), restart
  `coef_<name>_restart.nc` (with the RNG state: noisy coefficients continue identically).
  Coefficients are named by their symbol in the equation (`alpha`, `beta`, ...); the
  record a symbol loads is the setup's, in its par file. Terms are named by what they add.
  `tests/test_coefs.f90` covers the derived coefficients.
- `input/forcing_db.nml`: the library of named forcing groups that par files share:
  `snapshot_picontrol`, `_lgm`, `_b18_base`, `_b18`, `_b20`; `coef_alpha` (the default
  record, combined), `_dTo`, `_ramp`, `_hyst`, `_dT`, `_obm_north`, `_obm_south`;
  `term_lgm` (lgm - picontrol on alpha), `_b18`, `_b20`, `_ocean`, `_ramp`, `_hyst`,
  `_dT`, `_obm_north`, `_obm_south`; elsa groups `elsa_grl_lgp`, `elsa_grl_deglac`,
  `elsa_ant_paleo` (`[ytrc] elsa_nml`, `elsa_group`). Named groups are looked up in the
  par file first, then in the `var_defs` files in order (a par file that defines a group
  wins; a group defined nowhere stops the run, naming the group and the files).
  Experiments select with `-p forcing.coefs==... forcing.terms==...`.
- `input/clim_db.nml`: shared climate database of the snapshot backend, one varslice group
  per product field. Present-day references regional (`{domain}/{grid_name}`): MAR v3.14
  1961-1990, RACMO2.3p2 1981-2010, ERA5 1981-2010, ORAS5 1981-2010, ISMIP7 ocean.
  Paleoclimate from the global FesmData products, remapped online (`bilin`): PMIP3 and
  PMIP4 ensemble means, CLIMBER-3α, CLIMBER-X (with the anomalies to piControl),
  atmosphere and ocean; Buizert et al. (2018) `b18_*`, Badgeley et al. (2020) `b20_*`;
  PaleoMIST and GLAC-1D ice thickness. Path keywords `{grid_global}` (`[snapshots]
  grid_global`, default `GLOBAL-1DEG`) and `{exp}` (the snapshot's `exp`).
- `input/yelmox_defaults.nml`: defaults of the yelmox groups `[forcing]`, `[coef]`,
  `[term]`, `[snapshots]`, `[snapshot]`, `[taper]`, `[bias]`, `[orbit]`, `[ice_bnd]`,
  `[obm_forcing]` (templates of the named groups, domain suffixes); a par file needs only
  the keys that differ. Unknown keys stop the run (`nml_validate`).
- `.runme/info.json`: `par_defaults` (`input/yelmo_defaults.nml`,
  `input/yelmox_defaults.nml`), so `-p` can set a key of a group the par file lacks (e.g.
  `-p bias.active=True`), and `par_library` (`input/forcing_db.nml`): `-p` on a key of a
  library group copies the whole group into the staged par file, then applies the
  override. Needs runme >= 0.5.26.
- Climate transforms after the backend and the terms, in a fixed order, contracts checked
  at init, refused for esm with its smb and for rembo (`libs/climate_mods.f90`):
  `[taper<sfx>]` (the atmosphere towards the reference with a static weight field, `now =
  ref + w*(now - ref)`, the reference at the current surface; `outside = "mean"`: where `w
  < 1` towards the reference plus the mean anomaly of the cells with `w = 1`, per month),
  then `[bias<sfx>]` (additive temperature and multiplicative precipitation fields,
  monthly or annual, varslice groups with online remap; groups with a time axis are sliced
  at the forcing time, `bias_dtas_ann_now` in the 2D output). Both written to the climate
  2D file.
- Prescribed ice `[ice_bnd]` (`libs/ice_bnd.f90`, optional): Yelmo's `H_ice_ref` in its
  fixed and relaxation cells (e.g. the `yelmo_masks.border_width` strip) from an ice
  history, `method = "file"` (at the forcing time, clamped to the file's range) or
  `"scaled"` (a coefficient `coef` times one slice); `step_ice_bnd` before `step_icesheet`
  and at the cold start. `H_ice_ref` is in `yelmo.nc`.
- chion as surface model (`[comps] surface_method = "chion"`, `libs/surface_chion.f90`;
  `[surface_chion]` host keys, `[chion]` with `model = "itm" | "pdd" | "bessi"`): daily
  annual cycle on `grid_surface` over land and ice columns, host-supplied insolation, PDDs
  (sampled every 10 days) and snow fraction; tsrf and albedo are annual means of chion's
  surface state, the firn warming by refreezing from `chion_tsrf_annual`; the domain's
  physical constants; spin-up, restart, `chion.nc`. PDD: `pdd_method = "cg05"` (Calov and
  Greve, 2005), `[pdd_chion]`. BESSI: the host supplies the surface shortwave (TOA
  insolation x `[surface_chion] trans_sw`, and the TOA insolation for chion's cloud
  proxy), constant `wind_speed` and `rel_hum`, the air pressure from the surface
  elevation, the solar longitude and `H_ice`. `input/chion_defaults.nml`.
- Output precision: the 2D/3D fields of the yelmox writers (`isos.nc`, `mshlf.nc`, the
  surface file, the climate file, `yelmo_cmip.nc`) are written with a fixed precision per
  field (ncio `prec_abs`/`prec_rel`, bit rounding for compression): 1e-3 m elevations,
  1e-3 K temperatures, velocities 1e-3 m/yr or 1e-4 relative, 1e-4 m/yr mass balance (0.1
  mm w.e./yr), 1e-4 fractions, 1e-4 relative otherwise; masks, time series, restarts and
  static fields lossless. `[ctrl] write_lossless` (optional, default False): all output
  lossless (e.g. for bit-for-bit comparisons).
- `[comps] dt_isos` [yr]: isostasy update interval (10 in all par files; was every step).
  FastIsostasy integrates the interval in one call; the bedrock is held in between. Its
  yearly step was ~37% of the ANT-16KM ISMIP7 spin-up wall time. Answers change slightly.
- Timing of the time loop: the wall time of each part (output, forcing, relax, optimize,
  isostasy, coupling, ice sheet, climate, surface, shelf; bipolar also the ocean box
  model) printed at the end of a run and written to `timing.txt` at the `tm_1D` cadence.
- `yelmox_bipolar`: the ocean box model alone (`active_obm` without `active_north`/
  `active_south`); `atm2obm`, `ism2obm` and the ocean terms then stop the run.
- `[sim] recon_time`: the slice of `recon_path` for the `recon` cold start (was slice 1).
- `[smbpal] time_equil`: length of the smbpal ITM cold-start spin-up (was hard-coded 100
  yr).
- `{domain}`/`{grid_name}` in the isostasy `rheology_file`/`mask_file` (`grid_isos`), the
  `smb_simple` target mask (`grid_surface`) and the bipolar hydro masks (each domain's
  `grid_ice`); the Antarctic, LIS, bipolar and SRG par files use them (the ANT-32KM
  configs and bipolar south read the ANT-16KM Earth structure before).
- Greenland paleo on `yelmox/yelmox_GreenlandPaleo.nml` (GRL-PAL-16KM; docs page
  `experiment-grl-paleo.md`): ERA5 + ORAS5 present day with the MAR v3.14 - ERA5 1981-2010
  bias over Greenland (`[bias]`; `scripts/grl-paleo/make_bias_mar_era5.jl`), outflow-only
  "infinite" borders (`[yelmo] experiment = "infinite"`, `mask_border = "dynamic"`; needs
  the Yelmo `mb_border` budget, `input/yelmo-variables-ytopo.md`), prescribed border ice
  (`border_width` + `[ice_bnd]`, PaleoMIST) as an option, off; calving `vm-m16`, basal
  friction scaled with bed elevation, deep-ocean melt only in the open ocean below 2000 m,
  GHF Fox Maule (2005), isostasy LV-ELVA with two uniform mantle layers, sea level of
  Spratt and Lisiecki (2016) with the Holocene offset removed linearly over 0-10 ka
  (`input/sealevel_spratt2016.dat`, `make_sealevel_spratt2016.jl`). Setup A: spin-up and
  last glacial period on `alpha_combined_125kyr` (PMIP3 piControl/LGM, transient
  insolation, elsa); setup B: LGM spin-up from GLAC-1D and deglaciation after Tabone et
  al. (2024) (Buizert temperature, Badgeley precipitation, `[taper] outside = "mean"`
  beyond the GRL box, `taper_b18_b20` from `make_taper_recon.jl`). Scripts
  `scripts/grl-paleo/run_spinup.sh <forcing> <time>`, `run_lgp.sh`, `run_deglac.sh`. No
  optimization.
- Antarctic paleo (32 km) on `yelmox_Antarctica.nml` with `-p` (`scripts/ant-paleo/`): 15
  kyr optimization spin-up and -130 kyr to +10 kyr on the glacial index
  `alpha_combined_125kyr_interp`, sea level `sealevel_rohling_450kyr.dat` on relative
  time, ages traced with elsa.
- Experiment scripts on the ESM bases: `scripts/ismip7/` (optimization spin-ups
  `opt_ant.sh`, `opt_grl.sh`, `opt_grl_ismip.sh` by L. Gutierrez Gonzalez, and
  `opt_grl_k24.sh`, Greenland with the K24 basal hydrology: `method_transport = 1`,
  `k24_sliding_law = 4`; projections `run_ismip7_*.sh`; Antarctica on ANT-16KM),
  `scripts/1pctCO2/` (Greenland, Antarctica; the global CMIP fields of `[esm] esm_name`,
  default MPI-ESM1-2-LR, remapped online), `scripts/tipmip/` (Greenland; stabilisations
  branch off the ramp's restarts), `scripts/nudge/run_nudge_antarctica.sh`. Ice sheet and
  isostasy on; spin-ups pass `comps.dt_clim=10`.
- `scripts/run_chunks.sh`: a long run as a chain of short Slurm jobs (each chunk restarts
  from the bundle of the previous one; `afterok` dependencies; overrides for the first or
  later chunks); docs page `run-chunks.md`.
- Forcing builder (docs page `forcing-builder.md`, `scripts/forcing-builder/`): a base
  plus `-p` resolved, the forcing equation and its coefficients drawn over the run, the
  run written out as par file and runme command; Julia `export.jl`/`explain.jl`.
- `scripts/docs-figures/` (climate-only runs of the test cases, the forcing figures of the
  docs) and `scripts/climate-compare/` (climate-only comparison of two par files: domain
  means and maps of `t2m_ann`, `pr_ann`, `to_ann`).
- Docs: experiment pages (Greenland, Antarctica, Antarctic paleo, Greenland paleo, ISMIP7,
  1pctCO2, TIPMIP, Laurentide/North America, mountain glaciers), configuration groups,
  climate backends, run chunks; existing pages checked against the code.

### Changed
- Dependencies: fesm-utils, yelmo (with FastHydrology), FastIsostasy, rembo1, chion and
  insol at the v2.4-era versions pinned in `.configme/manifest.toml` (tags set at
  release). Existing builds need `make clean` once (stale `htopo.mod`,
  `insolation.mod`, `interp1D.mod` in `libyelmox/include`).
- One driver for the regional runs: `yelmox_esm` and `yelmox_rembo` are removed; `yelmox`
  chooses the climate backend at runtime, `[comps] climate = "snapshots" | "esm" |
  "rembo"` (replaces `make CLIMATE=`; `yelmox_snapesm.x` gone). REMBO is linked with `make
  yelmox rembo=1` (without it a stub stops the run). Runme `-e yelmox` (the `esm`, `rembo`
  aliases are gone); scripts in `scripts/ismip7/`, `tipmip/`, `1pctCO2/`, `rembo/`. Docs:
  `yelmox` and its bipolar variant instead of flavors (old pages redirected).
- Kryos modules: `libs/yelmox_domain.f90` split into `kryos` (domain type `kryos_domain`,
  was `ice_domain`; configuration, init, remap), `kryos_regions`, `kryos_coupling`,
  `kryos_startup`, `kryos_output`, `kryos_forcing`. The drivers write the coupling
  sequence out in the time loop (`yelmox_step` and the bipolar `advance_*` wrappers gone):
  `step_relax`, `step_optimize`, `step_isostasy`, `step_climate(force)`, `step_surface`,
  `step_shelf`, `couple_to_yelmo` (its own step before `step_icesheet`),
  `couple_yelmo_to_htopo`. The driver passes `dtt` to `step_optimize`, `step_climate` and
  `domain_startup`; the climate initializes at the forcing time `time_rel`.
- Time loop: output and restarts at the top of the loop for the current time, the final
  state and a final restart bundle on the last pass (also when a coefficient's kill switch
  trips). The zero-length first step at `time_init` is gone (it ran an extra optimization
  update and, in bipolar, an extra OBM step; cold starts change slightly). Restart
  cadence: `[tm_rst]` timeout group (`method = none|const|file|times`) replaces
  `dt_restart`. Timeline from fesm-utils `tstep_init`.
- Configuration groups: `[coupling]` is split into `[comps]` (`with_ice_sheet`,
  `with_isostasy`, `with_climate`, `with_surface`, `with_shelf`; `climate`,
  `surface_method`, `dt_clim`, `dt_isos`) and `[sim]` (cold start `init_method` = `none`,
  `equil`, `recon`, `recon_ref` with `init_equil_time`, `recon_*`; `init_marine_H`,
  `init_kill_shelves`, `init_time_thrm`; `relax`, `opt`, `use_negis`, `scale_glacial_smb`,
  `lim_pd_ice`); bipolar `[comps_<sfx>]`, `[sim_<sfx>]`. The cold start is the same in
  every driver (`domain_init_ice`; was chosen by domain name, and by `tstep_method` for
  Laurentide/North). `[ctrl]`: `run_step` (the group holding the timeline),
  `calendar`/`calendar_ref`, `restart` (the bundle to start from; bipolar one bundle per
  run, `restart-<kyr>-kyr/` with the shared bsl, obm and coefficient restarts and a
  subfolder per domain). Components named by the boundary they supply: `surface` (mass
  balance + temperature) and `shelf` (shelf-base melt + temperature): `smb_method` ->
  `surface_method`, `with_marine_shelf` -> `with_shelf`, `write_smb`/`write_mshlf` ->
  `write_surface`/`write_shelf`; `with_surface` and `with_shelf` need `with_climate`. The
  topography relaxation is its own switch `[sim] relax` with `[relax]` (`topo_rel`,
  `tau1`, `tau2`, `time1`, `time2`, `m`; was `[opt] rel_*`); after `time2` the `[ytopo]`
  values apply again. Old par files stop with "parameter not found".
- `[domain]` defines the domain (`[domain_north]`/`[domain_south]` in bipolar): `name`,
  `phys_const`, the grid of every component (`grid_hub`, `grid_ice`, `grid_isos`,
  `grid_clim`, `grid_surface`, `grid_shelf`; blank = default), `basins`. It replaces
  `[htopo]`, the `grid_*` keys of `[coupling]` and `domain`/`grid_name` of `[yelmo]`.
  Yelmo's grid comes from `maps/grid_<grid_ice>.txt`, and the hub topography (`[topo]`),
  remapped conservatively to `grid_ice`, is its initial topography and present-day
  reference (par files drop `yelmo.grid_path`, `yelmo_init_topo.init_topo_*`,
  `yelmo_data.pd_topo_*`). The hub reads the bed roughness `z_bed_sd`, fills the gaps of
  its topography file, keeps its reference geometry and adds Yelmo's bed displacement and
  ice-thickness change on a finer hub (`couple_yelmo_to_htopo`, was overwritten with
  Yelmo's fields); the hub geometry reaches `grid_clim`/`grid_surface` conservatively when
  finer; no Yelmo<->hub maps for identical grids. Answers change where `grid_ice /=
  grid_hub` (Antarctica, bipolar south: Yelmo 32 km, hub 16 km).
- Regions from FesmData v2 (`[topo]`, `[regions]`; fesm-utils `topodata`, `regions`; htopo
  moved to fesm-utils): each component makes its own masks from the regions on its grid
  with selection expressions: Yelmo in `[yelmo_masks]` (where ice is allowed, relaxation,
  error metrics, basins, named regions), the shelf in `[marine_shelf]` (`mask_ocean`,
  `mask_deep_ocean`, `mask_pico_deep`, `basins`, the corrections
  `corr_method`/`corr_names`/`corr_values`/`corr_<name>`, replacing `tf-grl`, `tf-ant`,
  `basin_number`, `basin_*_corr`), the climate in `[esm] basins` (ISMIP7 basins as a
  `[regions] layers` entry), the cold start in `[sim] recon_regions`, NEGIS in `[negis]
  region_*`. `libs/ice_sub_regions.f90` and the per-domain masks of the driver are gone.
  Antarctica: corrections `b1 b12 b14 b15` by Zwally basin.
- The hemisphere of a domain (seasons, lapse rates, esm forcing) follows from the latitude
  of `grid_clim`, no longer from the domain name.
- Physical constants: the domain owns them, `[domain] phys_const` (a group of
  `input/yelmo_phys_const.nml`, taken from yelmo: fesm-utils `phys_constants`, `rho_asth`,
  `rhoc_rock`), loaded once and passed to the hub, Yelmo, isostasy, the shelf, PICO,
  `smb_simple` and chion (`chion_phys_const.nml` and the density keys of
  `&isos`/`&smb_simple` dropped). `rho_uppermantle` (3400) is replaced by `rho_asth`
  (3300): equilibrium deflection ~3% deeper. The year is the CF/UDUNITS 31556926 s in
  every group (Earth was 365 d). The CMIP output conversions come from Yelmo's constants
  (mass rates -1.66 %, masses -0.76 %). Answers change.
- Snapshot backend `snapshots` (was snapesm; `libs/snapshots.f90`): reads its snapshots on
  the climate grid (online remap of `var_defs` groups, `remap`, `grid_src`); fields
  optional per snapshot, the reference sets what is built (atmosphere `tas`+`pr`+`zs`
  and/or ocean `to`, `so` optional: salinity 0 without); ocean annual or a monthly
  climatology; `tas_season = "jja" | "djf"` for two-field `tas`; Gaussian smoothing at
  loading `clim_sigma`, `ocn_sigma` [km] (PMIP3: 50 km); one monthly lapse rate for the
  sea-level reduction and the inflation to the surface, so a snapshot at its own elevation
  is recovered exactly (monthly `tas` changes by up to ~0.4 K per km of elevation in
  spring/autumn, summer `tas` and smb slightly). `snapshots_end`, `climate_end` in the
  drivers.
- The regional configs run on snapshots + terms (snapclim removed): Greenland
  (`pd`/`picontrol`/`lgm` on `alpha_orbital_r2020` normalised, ocean `term_ocean`),
  Antarctica (`pd`), North/LIS (PMIP3 on `alpha` = 1, `dTo` factor 0.25), Pyrenees (`pd` +
  a uniform term, `alpha_orbital_t2018` normalised 0 at present and 1 at 21 ka, `f_ta =
  -10` K), SRG (`pd`), bipolar (north as Greenland, south present day). Answers change:
  Antarctic reference RACMO2.3p2 1981-2010 t2m/pr (was the RACMO2.3/ERA-Interim hybrid:
  -1.0 K, +2.5% precipitation over the ice) and the ISMIP7 ocean (was ISMIP6 J20: +0.1 K);
  PMIP3 from the global product (GRL +0.3 K over the ice); the ocean anomaly from the
  sea-level air anomaly without lag (North/LIS LGM ocean 0.8/1.8 K warmer); dropped
  `f_hol` (GRL), the 80 km smoothing of the reference before 1 ka (`f_corr`), the SRG and
  Pyrenees oceans (`with_shelf = False`). Breaking for v2.3 par files: `[snap]` ->
  `[snapshots]` (`lapse`, `f_p`, `grid_global`; `ref_name` -> `[forcing] reference`,
  `var_defs` -> `[forcing] var_defs`; `combine`, `manifold`, `f_p_ne`, `f_stdev`, `f_hol`,
  `f_to`, `dTa/dTo/dSo_const`, `snapshots`, `indices` gone); `[snap_<name>]` ->
  `[snapshot_<name>]`; `[snap_idx_<name>]` and `[tsgen]` -> `[coef_<name>]` (the tsgen
  keys); `[snap_field_<name>]` -> the `snap0`/`snap1`/`fields` of a `[term_<name>]` (the
  rule is fixed by field); `blend = "fraction"` -> the derived coefficient `dTo` +
  `term_ocean`; `[tsforcing]` (`f_ta`, `f_to`, `f_so`) -> a `uniform` term; snapclim
  (`atm_type`, `ocn_type`, `[snap_hybrid]`, `[snap_clim*]`) -> snapshots on
  `input/clim_db.nml`.
- Climate on a coarser grid than the surface: `tas`/`pr` are moved from the climate grid's
  surface elevation to the surface grid's with the backend's monthly lapse rate and
  precipitation sensitivity. Unchanged when `grid_clim = grid_surface`.
- Climate output: each backend writes its own file, `[output] write_clim` (was
  `write_snap`): `snap.nc` (snapshots: the applied monthly `tas`/`pr`, `tsl_ann`,
  `z_srf`), `esm.nc` + `esm_ts.nc` (the forcing means on the climate grid; `t2m_1d`,
  `pr_1d` of the applied climate, `pr_1d` only when the climate has pr; was
  `yelmo_ts_esm.nc`), `rembo.nc` + `rembo_ts.nc` (`dT_sum`, `dT_ann`, `smb_mean`, `aar`).
  The 2D file has `t2m_sum` (summer mean of the applied `tas`) and the annual ocean
  profiles `to_ann`, `so_ann` for every backend that has them. CMIP/ISMIP output is a
  general option, `[output] write_cmip`/`dt_cmip` (was `[esm] write_formatted`/
  `dt_formatted`; `libs/cmip_output.f90`). The output switches need their component.
- ESM forcing is a climate backend (`climate = "esm"`, `[esm]` and its periods from the
  `run_step` group): the ocean at the shelf base (`T_shlf`/`S_shlf`, taken by the shelf in
  place of depth profiles; shelf base `z_srf - H_ice` as for the other climates, ESM
  answers change under the ice shelves), the surface mass balance (`surface_method =
  "climate"`: reference smb + anomaly, corrected from the present-day surface with the smb
  gradient; replaces `[esm] use_smb`) and subglacial discharge (`Qd`, read only if the sgd
  group's `filename` is not `"none"`). The ocean is extrapolated into the ice-shelf
  interiors on copies. `experiment = "ctrl"` in a transient runs on the reference;
  `esm_ant.nml`/`esm_nudge.nml` reference from the RACMO2.3p2 1979-2022 and ORAS5
  1979-2024 series (mean over `time_ref`). ISMIP7 runs `surface_method = "climate"` (the
  ESM smb of the protocol). `[esm] grid_src` and `&basins` removed.
- REMBO is a climate backend (`climate = "rembo"`, `surface_method = "climate"`): par
  files `yelmox/yelmox_rembo_Greenland.nml` and `yelmox/rembo_Greenland.nml` (staged by
  runme), script `scripts/rembo/`. REMBO supplies the atmosphere and smb; its summer
  anomaly `dT_sum` is the `tas` part of the uniform terms (the `pr` part of a term that
  also lists `pr`, e.g. `term_hyst`, is ignored with a notice at start-up; a uniform term
  on `pr` alone is refused); the ocean is an ocean-only snapshots reference (ORAS5; was
  none). The marine-ice rule and the cold start are those of `yelmox`.
- `yelmox_bipolar`: the OBM box air temperatures from one global anomaly coefficient,
  `[obm_forcing]` (`coef`, `gain`): `theta_i = theta_i_init + ampl_i × gain × f` for every
  box, vapour fluxes likewise (`libs/obm/atm2ocean.f90` gone). The ocean of each domain
  from the box model is the derived coefficients `obm_north`/`obm_south` and the `replace`
  terms (`[ctrl] obm2ism` gone); the box model and its coefficients run before the
  domains, so their startup sees the box ocean. The OBM restart is in the run bundle,
  keeps the `*_init` values of the par file and holds the box air temperatures and the
  overturning. `[forcing]` holds the coefficients, `[forcing_north]`/`[forcing_south]` the
  reference and terms, `[snapshots_north]`/`_south` the rules.
- Insolation from the [insol](https://github.com/fesmc/insol) library (`make
  insol-static`, `libinsol.a`; vendored `libs/insol/` removed). The orbit is
  `[orbit<sfx>]` (`libs/orbit_forcing.f90`), shared by chion's host and smbpal: `method =
  "forcing_time"` (at the forcing time), `"fixed_time"` (at `time`) or `"fixed_orbit"`
  (`ecc`, `obl`, `per`); it replaces `const_insol`, `const_kabp` and `insol_fldr` of
  `[surface_chion]` and `[smbpal]`. smbpal computes its insolation on the grid and only
  when the orbit changes. The insolation changes by up to ~2.5 W/m², mostly at low
  latitudes (insol integrates the daily mean exactly; the vendored copy summed 24 hourly
  samples).
- chion is the default surface model (`surface_method = "chion"` in the par files except
  LIS `smb_simple` and the ESM smb setups `climate`; the built-in default): ITM configs
  with `[itm_chion]` (`[chion] nml_itm`; Greenland `alb_ice = 0.31`, was 0.4, calibrated
  to MAR v3.11 1961-1990: runoff 635 Gt/yr, MAR 645), PDD configs (Antarctica, North,
  Pyrenees, bipolar south) with chion `pdd` (`cg05`) instead of smbpal; Greenland BESSI
  `[chion_const] alpha_ice = 0.55`. `[smbpal]` stays as an option. Answers change.
- Surface output: one file per surface model, `<surface_method>.nc` (`chion.nc`,
  `smbpal.nc` with smb in mm w.e./yr, `smb_simple.nc`); with `surface_method = "climate"`
  the climate file has the applied smb. The climate's smb is landed on Yelmo directly (no
  longer through smbpal); smbpal is initialized and restarted only when it is the surface
  model.
- smbpal ITM runs point by point with OpenMP, the daily forcing and insolation prepared
  once per call (ANT-8KM ISMIP7 10-yr test: 22.5 -> 8.0 min, 16 threads; equal to
  single-precision round-off). Daily output (`file_out_day`) removed.
- FesmData v2 inputs from the link `fesmdata` (`$FESMDATA`, copied into the run folders by
  runme): Greenland present day MAR v3.14 1961-1990 (reference climate, esm reference, pd
  smb/T_srf as skin temperature; was MARv3.11-ERA), North/LIS reference
  `Atmosphere-ERA5-1981-2010` (was ERA-Interim), ocean `Ocean-ORAS5-1981-2010` (was ORAS4;
  esm Greenland from GRL-16KM remapped), paleo `Paleoclimate-PMIP3-*` (was
  `PMIP3_sig50km`), velocities `Velocity-Joughin2018`/`-Rignot2011`, GHF
  (`GHF-Shapiro2004`, `-FoxMaule2005`, `-Loesing2021`, `-HazzardRichards2024`), sediments
  `Sediment-CRUST1` (was `SED-L97`; SRG keeps it), Antarctic basal melt
  `BasalMelt-Rignot2013` (`bmelt_eq_reese`, gaps filled from the nearest value),
  isochrones `Stratigraphy-MacGregor2025`, ISMIP7 Antarctica pd smb/T_srf from RACMO2.3p2
  1981-2010 and the esm reference from `SMB-RACMO2.3p2-1985-2014`, Greenland
  `TOPO-M17-v5`, SRG `TOPO-GLO90-GEBCO2026`, rembo1 on `TOPO-BedMachine-v6` and MAR v3.14.
  Pyrenees: ice sheet and surface on `PYR-250M`, climate (ERA5, was WorldClim), isostasy
  and GHF on `PYR-2KM`, `dt_clim = 50`, chion `time_equil = 10`. GHF and sediment obs are
  remapped onto the Yelmo grid from the grid of the file. GHF units from the variable's
  `units` attribute (W/m2 or mW/m2), missing values filled. Answers change.
- Follows yelmo dev (`input/yelmo_defaults.nml`, `yelmo_phys_const.nml` and the
  `yelmo-variables-*.md` tables synced): renumbered `ytherm.qb_method` (par files keep 2,
  now "faces to quadrature nodes"), capacity basal BC, `ytherm.method` the default
  `"enth"` (the par files set `"temp"`), `frz_scale`/`frz_efold`/`frz_min` (replace
  `slide_T`, `gamma_T`, `lambda_min`), `strain_heating`, `gl_temperate`, `pc_rho_max`,
  `pc_tol = 1`, `pc_eps = 0.02`, `de_max = 100`, `bkt_floating_mode = 0`, velocity-limit
  drag (`ssa_vel_lim_method`, `ssa_vel_max = 10000`), subgrid fronts (`front_subgrid =
  "marine"`, `tau_ice_flt`/`tau_ice_grnd`; `f_ice_method` gone), `gl_sep = 3` (ESM
  Antarctica, ISMIP7/1pctCO2 scripts; `gl_sep = 2` removed), `ssa_lat_bc = "all"` (SRG and
  ANT paleo were `"floating"`), `grad_lim_zb = 0.5` (was 0.05: ANT-16KM grounded thickness
  RMSE 131 -> 16-24 m at the start), SSA energy solve CG + SSOR, tol 1e-3, `tau_ice_flt =
  100e3`, `vm-m16-old`, `mask_border`, K24 keys renamed as in FastHydrology (e.g.
  `k24_N_ub_coupled`, `k24_water_viscosity`, `k24_flux_solver = 3`), `[ytrc]` elsa/tracer
  keys. Removed keys dropped from the par files (`dHdt_dyn_lim`, `surf_gl_method`,
  `margin2nd`, `ssa_beta_max`, `scale_T`, `T_frz`, `use_strain_sia`, `[opt]
  cf_min`/`opt_cf_min`/`use_yelmo_cf_min`). Answers change.
- `[opt]`: `opt_cf = "none" | "L21"`, `opt_tf = "none" | "L21" | "L21-points"` (were
  logical; `"L21"` one `tf_corr` per basin, `tf_basins`; `"L21-points"` each point, the
  former method). `step_optimize` only optimizes `cb_ref` and `tf_corr`; `tf_corr` is
  optimized on `grid_shelf` in place (Yelmo-grid inputs remapped there; the basins, and
  the `tf_basins` ids, are `[marine_shelf] basins` on the shelf grid), no longer through a
  con/bilin round trip through the Yelmo grid each step: unchanged when `grid_shelf =
  grid_ice`, answers change for Antarctica (shelf 16 km, ice 32 km). The optimization's
  cold-start `cb_ref` is set before `yelmo_init_state` in every driver (`opt` cold starts
  change slightly).
- NEGIS (`[sim] use_negis`) and the glacial smb scaling (`scale_glacial_smb`, group
  `[glacial_smb]`: `dt_lgm`, `lat_lim`, `fac_lim`) apply to any domain. NEGIS takes its
  parts from the `NEGIS` basin set (`[negis] region_centre`, `_south`, `_north`), sets
  `ytill.method = -1` and is refused with `opt_cf`; the `cf_0` -> `cf_1` ramp (-11 to 0
  ka) scales all parts; the factors depend on the forcing time only and no longer
  overwrite the parameters. Greenland no longer sets `cb_ref = cf_ref` for `ytill.method =
  -1` (a cold start without optimization starts from Yelmo's `cb_ref = 1`).
- Isostasy: the bedrock anomaly is coupled to Yelmo (`z_bed = z_bed_ref + w + we`; was the
  absolute bedrock on `grid_isos`, which smoothed every trough and sill finer than
  `grid_isos`), the ice load and reference coarsened conservatively (were bilinear), cold
  starts of every driver through `domain_init_isostasy`. `check_isostasy_reference` exits
  with status 1. Answers change.
- Model time in double precision (`real(dp)`): the drivers' `dtt`, `[comps]
  dt_clim`/`dt_isos`, `[output] dt_cmip` and the time arguments of the kryos modules, the
  climate, surface and shelf components and the OBM coupling.
- `scripts/check_par_nml.py` checks the par files against `input/yelmo_defaults.nml` and
  `input/yelmox_defaults.nml` (unknown keys, duplicate groups and keys, named groups
  against their templates, the snapshot bindings against the `var_defs` files, the
  library); `--dead` reports keys no source reads. Checks `yelmox/*.nml` and
  `yelmox_bipolar/*.nml` by default.
- Par files: `yelmox/` holds one base per domain and climate backend (Greenland,
  GreenlandPaleo, Antarctica, North, LIS, Pyrenees, SRG, esm_Greenland, esm_Antarctica,
  rembo_Greenland) plus the REMBO file `rembo_Greenland.nml`; an experiment is a `-p` list
  in its scripts (`scripts/<exp>/params.sh` when shared), run control in the invocation.
  Groups in sections following the time loop; production run lengths; keys and groups that
  are not read removed (e.g. `[ice_enh]`, unused `corr_*`, `mask_pico_deep` outside pico);
  no personal paths. Experiments that were par files are one-line commands in the docs
  (present day, ramp).
- Build: OpenMP by default (`openmp ?= 1`; `make <target> openmp=0` serial); yelmox
  objects depend on the archives of the dependencies (`libfesmutils.a`, `libyelmo.a`,
  `libisostasy.a`, `librembo.a`), so a type change in a dependency rebuilds them.
  Regenerate the Makefile with configme.
- Run scripts take the queue options from the environment if set (`runopts`,
  `HPCOPT_SPINUP`/`HPCOPT_SCEN`, `RUNOPT`); runme's `--mem` follows the config (was
  hard-coded 30Gb, which broke OpenMP jobs on Levante's shared partition).
- Yelmo's internal output (regional 1D files) goes to the domain's output folder; the
  numerics/speed metrics are in Yelmo's `yelmo_metrics.nc`.
- FastIsostasy `[isos] dt_method = "rkc"` (Runge-Kutta-Chebyshev) available;
  the par files keep `tsit54`.

### Fixed
- Exact restarts (bit-for-bit with the continuous run, `yelmox` and `yelmox_bipolar`): the
  isostasy restart carries the ODE solver's time step and state, the sea level is restored
  as saved (`bsl_startup` without `bsl_update`; each domain keeps its own restored bsl),
  the chion and smbpal restart files carry the annual fields of the last cycle (`ann_*`)
  and `domain_startup` continues the restored state, re-deriving only the climate
  (whatever `dt_clim`; a restart with `dt_clim <= 0` left it unset) and the stateless
  `smb_simple` surface (was a re-step of surface and shelf, one extra snowpack year).
  Bundles without `ann_*` restart inexactly, with a warning. A bundle restart initializes
  Yelmo's tracer backends (elsa, tracer; crashed in the first step). A coefficient without
  a restart file in the bundle keeps its cold-start state.
- `smb_simple` got the surface grid's coordinates in km labelled as metres, so its signed
  distance was 1000x too small; now in metres. LIS answers change.
- `ism2obm` freshwater flux (`calc_fwf`): ice volume converted to Sv with rho_w/rho_ice
  (21% too high), `[ctrl] fwf_definition` ignored (always `dVdt`), `dVdt` summed over the
  whole domain; now `dVdt`, `mb` or `ablation` within the hydrographic mask, in water
  volume; Yelmo fields no longer modified. Bipolar answers change with `ism2obm`.
- Marine shelf: freezing point at the ice-base depth below the current sea level (`z_base
  - z_sl`; was relative to present-day sea level), answers change when `z_sl /= 0`
  (paleo). `*-slope` methods: `quad-nl-slope` uses the slope-scaled gamma, base slope from
  the gradient magnitude (opposite x/y slopes cancelled). `bmb_method`, `tf_method` and
  their combination checked at init.
- PICO: box k fed by the basin-mean T/S of box k-1 and the basin box-1 overturning (Reese
  et al., 2018); box-n salinity sign; parameters from `[pico<sfx>]`; lakes excluded;
  `bmb_shlf` set in every cell (the basin's box-1 melt in partially grounded cells, for
  Yelmo's `bmb_gl_method`; 0 where fully grounded, on land and in lakes); no NaN without
  shelf beyond the grounding line. No shipped setup uses PICO.
- CMIP output: `hfgeoubed` in W m-2 (was x1e6 too large), `ligroundf` per area (kg m-2
  s-1); 2D fluxes from Yelmo's applied smb/bmb/cmb and the 1D integrals from them
  (`tendacabf` the applied smb, `libmassbfgr`/`libmassbffl` the applied bmb).
- ESM forcing: transient ocean anomalies were NaN/garbage for depth-less ocean forcing
  (Greenland ISMIP7, 2D monthly; fesm-utils `varslice_sub_mean`); `dT_shlf`/`dS_shlf` are
  the anomalies `dto + dto_var`/`dso + dso_var` for every domain (Greenland had the
  absolute temperature plus `dto`, `tf_shlf` ~274); the historical direct-SMB anomaly
  against the ESM's own reference-period mean; the reference ocean extrapolated on its own
  depth axis; a transient ocean field whose layout differs from its reference stops.
  Output: `pr_ann`/`pr_1d` the annual mean of `pr*dpr` in mm/d, `t2m_ann`/`t2m_sum` and
  the 1D anomalies annual and hemisphere-summer means (were January only), `dpr_ann` units
  `1`. ESM spin-ups need re-running.
- ESM par files: Greenland had the Antarctic ITM parameters (SMB strongly negative);
  ISMIP7 Antarctica spin-up `time_ref` outside the reference axis (missing climate); the
  ISMIP7 GHF (W/m2) was never converted and clamped to 0.1 mW/m2 (`obs_err_name` the GHF
  itself); `&ghf` keys and current topography/GHF files for esm Antarctica; meltMIP ocean
  reference from `ice_data/ISMIP7/`. ISMIP7 spin-ups must be redone.
- `opt.cf_init <= 0` starts the optimization from the till friction of the bed (`cb_ref =
  cb_tgt`) in every driver (was a negative friction except in rembo).
- `check_isostasy_reference` compares the reference bedrocks on the isostasy grid (the
  same bedrock agrees to round-off on any grid; SRG on 16 km failed the 10 m tolerance).
- `sed.method` is read (was uninitialised, so sediments were off).
- Geothermal: `ghf_err` filled and converted like `ghf`; the `ghf_min` floor applied
  whether or not the error field is read; the error field read only when `f_stdev /= 0`.
- smbpal PDD (`abl_method = "pdd"`): snowfall from `calc_snowfrac` on the monthly fields
  (was all precipitation). smbpal state zeroed at allocation (the PDD restart held
  uninitialized `alb_s`).
- Cold start: the surface is updated once more on the final initial geometry (new margin
  ice had smb = 0 until the next update); Yelmo topo diagnostics reconciled after the
  isostasy init.
- `[sim] scale_glacial_smb = True` is refused at start-up with `climate = "rembo"` (the
  glacial index needs the climate's reference annual air temperature, which REMBO does not
  supply; the run used an unallocated field).
- `[domain]`: `dx_hub` after `yelmo_init_grid` (was from the unset Yelmo dx).
- Output: `write_clim` without a climate crashed; marine_shelf restart field units (were
  all m/yr) and real `xc`/`yc` axes in km; OBM air temperatures in degC (were K), restart
  `m` in m3/yr (was Sv). nautilus `ft` initialised from `ft_init` (output only).
- rembo: `dT_ann` is REMBO's annual anomaly, (1+T_wintfac)/2 * dT_sum (needs rembo1
  `rembo_dT_ann`), instead of 1.3 * dT_sum (diagnostic).
- smbpal, itm, sediments and geothermal par loaders honour the value of `init`.
- SRG (Patagonia) runs again: `maps/grid_SRG-250M.txt` is the FesmData v2 grid (UTM zone
  18S, 250 m, 208 x 120); every component on `SRG-250M`; isostasy off by default (on:
  `SRG-16KM`).
- `maps/grid_GRL-32KM.txt`: 53 columns, as the data (the ncio size check stopped bipolar
  on the GRL-32KM climate). Delete cached `maps/scrip-*GRL-32KM*.nc` maps.
- `.runme/info.json`: the REMBO namelist is staged for `yelmox` runs.
- `tests/test_timestepping.f90`: current fesm-utils API, with checks.

### Removed
- snapclim (`libs/snapclim.f90`, its tests, `input/greenland_snp.nml`,
  `input/greenland_clim.nml`, `yelmox_Greenland_snapclim.nml`) with its `f_p_ne` and basin
  masks; `[tsforcing]`, `[tsgen]` (see Changed: the terms).
- Variant par files, replaced by `-p` lists on the bases: `yelmox_Greenland_snapesm`,
  `_snapesm_ramp`, `yelmox_pd_Greenland`, `yelmox_pd_Antarctica` (one-line commands in the
  experiment pages); `yelmox_esm_*` ISMIP7, K24, 1pctCO2, TIPMIP and nudge variants
  (arrays in `scripts/ismip7/`, `1pctCO2/`, `tipmip/`, `nudge/`).
- Legacy single-grid programs (`<flavor>/legacy/`, `make <flavor>-legacy`) and retired
  flavors (`yelmox_ismip6`, `yelmox_nahosmip`, `yelmox_rtip`) with their par files, make
  targets, runme aliases and `scripts/ismip6-2300.md`; the legacy ISMIP7 and ant-paleo
  scripts.
- Unused modules: `libs/ismip6.f90`, `libs/hyster.f90` (and the `yelmox_hysteresis_help`
  rule), `libs/simpleclim.f90`, `libs/gaussian_filter.f90`, the basal hydrology modules,
  `libs/stommel.f90` and the OBM scratch programs (`[ctrl] obm_name` must be `nautilus`),
  `timeline_init`.
- `[marine_shelf] tf_method = 0` (no external thermal forcing exists).
- Dead code in smbpal (`update_2temp`, `write`/`write_init`, `smbpal_end`, ...),
  sediments, geothermal, `surface_chion`, marine_shelf/PICO (c_deep smoothing, unused
  routines), the snapshots backend, `esm_forcing` (`esm_write_init`, `use_ref_atm/ocn`);
  `[itm] trans_c` (unused).

## [v2.3.1] - 2026-07-25

### Fixed
- `yelmox_rembo`: writing `rembo_ts.nc` aborted in ncio (`nf90_inq_dimid` for
  `xc`) because `yelmo_write_step_model_metrics` also emits the 2D `pc_tau_max`
  field, which needs `xc`/`yc` dims absent from the 1D timeseries file. The
  metrics call is disabled for now; numerics/speed diagnostics will move to a
  dedicated `yelmo_metrics.nc`.

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
