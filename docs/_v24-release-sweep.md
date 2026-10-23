# v2.4 release sweep (2026-10-23)

Status: implemented on branch v24-sweep: A-D, E33-E34, H, E25, E27, E28 (E29
left by decision); the release surface (A, H docs) on v24-release. F and G remain
post-release.

Read-only review of yelmox dev 809e179 + yelmo dev 309ff11c by ten parallel Opus
agents (nine code areas, one release surface). Nothing edited. The headline
claims below were re-verified by hand in the source. File:line refer to the
state above. `yelmo/` = `/Users/alrobi001/models/yelmo/src`.

## Verdict

The benchmarked default paths (chion surface; anom/lin/quad-nl shelf melt;
DIVA + energy solver + beta_method=1; enthalpy; impl-lis advection; bucket
hydrology) show no bug that changes answers. The mass budget closes and every
term is booked once. The CHANGELOG claims hold apart from a few stale lines.

Not ready to tag as is. Four release-surface items block, three real bugs
change answers in shipped configurations (LIS, SRG/PYR, bipolar), and two
documented options cannot start. All are small fixes.

## A. Release blockers (release surface)

1. **`input/yelmo_defaults.nml` out of sync with yelmo dev.** yelmox copy has
   `tau_ice_flt = 250e3` and lacks `vm-m16-old`; yelmo has 100e3. Runs read the
   yelmox copy. Shipped par files all set tau_ice_flt, so only external par
   files are hit. CHANGELOG:431 "synced" claim is false. Re-sync.
2. **configme manifests pin `2.4-beta`/`dev`** (`.configme/manifest.toml`):
   2.4-beta is behind dev in every repo (yelmox 195, yelmo 18, fesm-utils 16,
   FastIsostasy 4, rembo1 3, chion 3 commits); dev calls `isos_init(domain=,
   grid_name=)` which 2.4-beta FastIsostasy lacks. Tag each dependency and pin.
3. **yelmox main has 4 commits not on dev**: 803bf8d (v2.3.1 rembo fix +
   CHANGELOG section), 04eaa39 (experimental `YELMOX_TF_BASIN` tf_corr), two
   gitignore commits. Files they touch no longer exist on dev. Decide: add the
   v2.3.1 CHANGELOG section to dev, drop 04eaa39.
4. **yelmox CHANGELOG [Unreleased] narrates intermediate states**: `### Added`
   twice (8, 434), `### Changed` twice (164, 482); lines 49-62 describe
   `snap_*`/`idx_*`/`[indices]` (now `snapshot_*`/`coef_*`/`term_*`); 222
   `[idx_<name>]`; 229-233 old orbit method names; 643 `[coupling] climate =
   "snapclim"`; 427-428 snapclim dp changes. Consolidate to final state.
   yelmo CHANGELOG: lines 5-166 are uncategorised bullets above the `###`
   sections; line 43 says border reset booked in `mb_resid` (code and line 18:
   `mb_border`); ~1165 says redo restores only tpo/dyn (also hyd); 485 stale
   `ssa_lat_bc` default; two `yelmo_data` default entries disagree (54 vs 73).

## B. Bugs that change answers in shipped configurations

5. **smb_simple gets km coordinates labelled "m"** (`libs/kryos.f90:374-376`,
   `units="m"` with `grid_s%x` in km). Signed distance 1000x too small: inland
   accumulation decay off, plastic profile collapses. LIS runs since the
   multigrid port (5f6dbf2). Pass `grid_s%x*grid_s%cs%xy_conv`. Then check
   whether LIS `c0`/`L_acc` were tuned with the bug.
6. **`yelmo_init_grid_fromgrid` drops `k0`** (`yelmo/yelmo_grid.f90:363-367`):
   coords defaults k0=1; SRG/PYR (transverse Mercator, k0=0.9996) get lon/lat
   and `grd%area` ~0.08% off, and the written CRS says k0=1. Pass `k0` (or use
   coords `grid_init_from_grid`).
7. **Bipolar freshwater flux ignores `fwf_definition` and the hydro masks**
   (`libs/obm/ice2ocean.f90:66` hard-codes `fwf_def = "dVdt"`; masks read in
   `obm_coupling.f90:159-170` only used by the unreachable `dVdt_mask`). The
   `ablation` branch also clamps the live tpo arrays in place (intent-less
   allocatable dummies). Make `calc_fwf` pure, honour fwf_def, mask every branch.

## C. Documented options that stop or crash at init (no shipped config uses them)

8. **`[ice_bnd] active = True` stops**: `ice_bnd.f90:71` reads `var_defs` with
   `defaults_group="ice_bnd"`, key removed from defaults by 2035098. Use
   `forcing_var_defs(filename)` like taper/bias.
9. **Restart with `dt_clim <= 0` never computes the climate** (`kryos_startup.f90:98`,
   `kryos_coupling.f90:404`): `clim%now` stays unallocated. Add a `force` path
   separate from cold-start init.
10. **`write_clim` with `with_climate = False`** writes unallocated `ta_ann`
    (`kryos_output.f90:119`); **coefs with `write_yelmo = False`** append to a
    yelmo_ts.nc that was never created (`:75`). Resolve effective switches once
    in `domain_ctl_load`.
11. **`[bias]` + esm with its own smb**: `bias_apply` sums unallocated `pr`
    (`climate_mods.f90:436`). Extend `climate_caps_class` and refuse at init.
12. **NEGIS does nothing if switched on** (`kryos_coupling.f90:369`,
    `kryos_regions.f90:61-92`): cb_ref overwritten by `calc_ydyn` unless
    `till_method = -1`; `cf_0/cf_1` ramp computed, never applied; `cf_south/north`
    mutated permanently. Make it a proper external-cb_ref mode or remove.

## D. Wrong output / diagnostics (no effect on the model state)

13. CMIP `hfgeoubed` x1e3 should be x1e-3 (`cmip_output.f90:295`; Q_geo is mW m-2);
    2D `ligroundf` missing `/dx` (`:153`); 2D vs 1D use different smb/bmb fields.
    ISMIP7 submission variables.
14. `esm_ts.nc pr_1d` uninitialised when esm supplies the smb
    (`yelmox_climate.f90:594`; ISMIP7 GRL).
15. rembo `dT_ann = 1.3*dT_sum` stale: `T_wintfac = 2.0` gives 1.5
    (`yelmox_climate.f90:853,872`). Take it from REMBO.
16. `term_<name>` output is 0 for rembo (dup term state in `snap_term_class`,
    `snapshots.f90:158`); rembo snapshot terms never written.
17. Missing observed velocities enter rmse as 0 not mv (`yelmo/yelmo_data.f90:346`);
    `pd_err_z_bed` never computed (`:610`).
18. Applied-rate fields (dHidt, mb_net, mb_resid, cmb, ...) stale when
    `topo_fixed`/dt=0, and the restart reads them (`yelmo/yelmo_topography.f90:556-651`).
    Same class as the `mb_relax` fix already in the CHANGELOG. `adv_lin_fail`
    likewise (`yelmo_ice.f90:520`).
19. Always-zero fields in tables/C API/restart: `mat%now%strs`, `dQsdT`,
    `dyn%now%de_eff`, `duxdz_bar/duydz_bar`, `tpo d{zs,Hi,zb}d{x,y}_aa`;
    `ssa_err_acx/acy` units "m/yr" should be "1".
20. marine_shelf.nc labels T/S/tf fields "m/yr"; xc/yc axes 0..nx-1 "km"
    (`kryos.f90:405` passes no axes); OBM T "K" but degC; restart `m` "Sv" but m3/yr.
21. `nautilus%ft` initialised from `tt_init` (`libs/obm/nautilus.f90:114`); ft unused.

## E. Inconsistencies in non-default options (fix or mark unsupported)

22. `ytherm.dt_method` AB/SAM blends the already-blended field, not AB2
    (`yelmo/physics/thermodynamics.f90:528,603,...`). Remove or keep raw rates.
23. `t_dep_source="trc"` isochrones treat MV=-9999 as a deposition time
    (`ice_tracer.f90:813`).
24. `method="temp"` drops hydrology heat `Q_wat` that the enthalpy path includes
    (`ice_enthalpy.f90:160 vs 774`); the basal-BC block is duplicated.
25. Shelf freezing point uses `z_base` (datum) while T/S are interpolated at
    depth below current sea level (`marine_shelf.f90:256-261` vs `513`); ~0.09 K
    at LGM. Convention decision: pass `z_base - z_sl`.
26. marine_shelf restart overwrites par-derived `bmb_ref`, `bmb_corr`,
    `tf_corr_basin`, `mask_ocn_ref` (`:863-875, 2038-2055`) while basins are
    recomputed from the par. Restart only true state.
27. `tf_corr` makes a con->bilin round trip every opt step when shelf grid !=
    ice grid (`kryos_coupling.f90:89-104`; ANT 16 vs 32 km): smoothed every step.
    Keep it on the Yelmo grid. Changes ANT opt answers.
28. K24 hydrology state (`k24_routed`, `k24_abs_g`, `k24_phi0`, `diss_face`) not
    in the restart (`yelmo_io.f90:1063`); restarted `method_transport=1` runs
    not exact. The restart reader is a hand-written second list. Either add the
    fields or state the limitation in the CHANGELOG.
29. pc mask isolated-point test compares raw `pc_tau` [m/yr] with `pc_eps`
    (scaled-eta target; `yelmo_timesteps.f90:288`): the 1.0->0.02 change lowered
    this threshold 50x. `pc_dt(1)` stores clamped dt (`yelmo_ice.f90:490`).
    Both answer-changing; defer until a re-benchmark.
30. `tf_method=0`, `*-slope` melt (gamma2D not passed; signed slope sum),
    PICO (wrong par group, lakes, stale melt, 0/0 into NINT, not Reese 2018),
    stommel OBM (tendencies /dt): unused and broken. Mark experimental or remove.
31. `taud_gl_method` unvalidated, `case DEFAULT` silent (`velocity_general.f90:1348`);
    `zeta_scale` DEFAULT silently linear; `calc_zeta` comment lists "wave".
32. `ssa_lis_opt_residual` default is minres for a non-symmetric matrix.
33. `strn2D%f_shear` allocated, never deallocated or zeroed (`yelmo_dynamics.f90:1027`):
    `yelmo_end` + `yelmo_init` fails. One line each.
34. `bnd%restart_relax_init` read every step, set only on restart (`yelmo_ice.f90:1489`).
35. Model time rounded to sp on the level-set path (`yelmo_topography.f90:295,935`;
    `lsf_module.f90:283`), contradicting the dp-time guarantee.

## F. Performance (no answer change unless noted)

36. **Redo snapshot deep-copies tpo/dyn/hyd every step** (`yelmo_ice.f90:173`),
    incl. dead rk4 arrays: ~0.4 s/step at ANT-8 measured standalone. Snapshot
    only when `dt_now > dt_min`; drop rk4.
37. **`pc_tau_max` running stats every step**, strided over 50 levels
    (`yelmo_ice.f90:525`): ~0.2 s/step at ANT-8. Only with `write_metrics`.
38. **Picard loop invariants recomputed per iteration** (DIVA default): active
    faces (twice), ATT^(-1/n), velocity-independent beta for `beta_method=1`,
    F-integral weights; diagnostics inside the loop. Energy assembler rebuilds
    sparsity and calls `assemble_row` twice per iteration
    (`solver_ssa_ac_energy.f90:187-319`). LIS init/finalize per solve
    (`solver_linear.F90:246`); `linear_solver_init` per advection solve
    (`solver_advection.f90:108`).
39. `update_ice_fraction` runs ~8x per pc stage (front cells, lsf area fraction,
    7-8 allocs each); `f_grnd_bmb` recomputed every stage.
40. Heap allocations inside per-column OpenMP loops (`ice_enthalpy.f90:651,1019`);
    `calc_dzeta_terms` recomputed per column in tracers (dup of the enthalpy copy);
    tracer column loop and `calc_magnitude_from_staggered` without OpenMP;
    `apply_tendency` serial.
41. Wrong loop order: `calc_advec_horizontal_column` k innermost
    (`thermodynamics.f90:231`); SSA visc loop i-outer (`velocity_ssa.f90:498`);
    SIA `tau_xz` 3D arrays unnecessary; `calc_pc_eta` j-innermost (roundoff).
42. yelmox: `smb_simple_set_mask` nc_read + O(N*Nb) distance every step
    (`kryos_coupling.f90:463`); static `H_ice_ref`/`z_srf` remapped every step
    (`:548,:417`); 3D ocean remapped clim->shelf every step instead of per
    dt_clim (`:594`); `find_open_ocean_and_lakes` sweep per update (capped
    qmax=1000, silent); basin means O(n_basin*N); time-varying snapshots
    re-reduce everything and copy `ref` every update; `chion_get_surface`
    fetches unused tsrf/albedo daily; `log_timestep` `nc_time_index` O(n^2).

## G. Structural cleanup (post-release candidates, answer-neutral)

43. SSA and DIVA Picard drivers are near-duplicate and have drifted (tau_b
    underflow 1e-5 vs 1e-15; `set_inactive_margins`; `visc_method=0`). One
    driver with a shear switch.
44. Boundary-string vocabulary parsed in >5 places with dead aliases
    (`EISMINT`, `periodic-xy`, `fixed`; EISMINT branch off by one); BC table
    copied in both assemblers and basal_dragging. Parse once to the code.
45. Present-day topography cleaned twice with different rules
    (`yelmo_data.f90:247` vs `yelmo_ice.f90:1197`): fixed/relax targets differ
    from the initial ice in thin cells. One routine (slightly answer-changing).
46. Calving: `eps_eff` from `dyn%strn2D` vs `tau_vm` from `mat%strn2D`
    (`yelmo_topography.f90:693 vs 704`); exp5 law implemented twice; aa->ac
    stagger block copied 3x.
47. Transport sums to zero only for impl-lis; expl/sico borders unbooked;
    `"zero"` borders drop inflow silently (`solver_advection.f90:306-421`).
48. Climate terms contracts: derived coef on a snapshots term lags one update;
    no single producer for derived coefs (bipolar both domains set `dTo`);
    replace discards snapshot terms on the field; unbound term fields ignored.
    All init-time checks.
49. Duplicated helpers: `lapse_monthly` (esm vs snapshots), `calc_snowfrac`,
    `solar_longitude_deg`, insol table (chion vs smbpal), `gauss_values` and
    magnitude in ice_optimization, cdf/error_function, per-lib precision kinds
    (`smbpal_precision`, sediments, geothermal, smb_simple, esm_forcing).
50. Heat capacity two sources (`thrm%now%cp` vs recomputed `cp_eff`); `Q_strn`
    from `mat%now%visc` not the dynamics viscosity (design decision); tracers
    use `H_ice` not `H_ice_dyn`; active-rock mask differs from ice column.
51. Restart: run-level objects (bsl, coefs, OBM) restored asymmetrically between
    drivers; `counter_restart` not restored (`kryos_forcing.f90:343`); restart
    reads thrown away right after (mask_bed, basins, regions, pc_tau*).

## H. Dead code and stale text (one housekeeping commit)

- Delete: `libs/stommel.f90`, `libs/gaussian_filter.f90` (buggy copy of
  fesm-utils), `libs/obm/zirajas.f90`, `libs/obm/run_nautilus.f90`,
  `yelmo/physics/solver_advection_new.f90` (compiled, unused),
  `yelmo/yelmo_regridding.f90` (empty stub), `yelmo/libs/basal_hydro_simple.f90`,
  all of `runge_kutta.f90` + `tpo%rk4` (`use_rk4 = .FALSE.`), `calc_l14_scalar`
  (precedence bug), `file_diagnostics` arg, `yelmo_get_precision` (wrong),
  `if (.FALSE.)`/`if (.TRUE.)` blocks in yelmo_ice/timesteps/velocity_*/
  deformation (one with a sign error), `ij = maxloc` per step.
- Dead parameters: `ytherm.solver_advec` (documented, no effect),
  `ymat.flow_law`, `ydyn.beta_gl_sep`, `use_ssa`, `[snapshots] combine` and
  `f_p_ne` (no-op; CHANGELOG:414 says removed), `ctl.tstep_method`, `ctl.dtt`
  dup, `[glacial_smb]` with no template, 14 unused `corr_*` in
  `yelmox_Antarctica.nml:718`, `mask_pico_deep` set in non-pico par files,
  smb_itm `trans_c`, `sed_limit`, smbpal `*_write*`/`_2temp`, `init_pars =
  present(init)` pattern (ignores the value) in 5 modules.
- Makefiles: `latinhypercube.o` rule (no source), `control.o` rule, commented
  rules, `yelmo_thermo` list with `interp1D.o`, obm deps missing `obm_defs.o`,
  `cleanall` skips chion, `usage` lacks test targets, `check_githash` unused.
- Docs: `README.md:7-13` pins fesm-utils 3f415cc and mentions ismip6;
  `configme-install-details.md` stale (4 repos, `fesm-utils/utils`);
  `submodules-plan.md` published and outdated; `database-namelists.md:10,239`
  ismip6; public pages link to unpublished `_forcing-terms.md`; `_forcing-design`
  and `_snapclim-retirement-plan` status lines say "not started";
  `tools/yelmo-config` defaults snapshot stale.
- Wrong comments (selection): qb_method numbers (`thermodynamics.f90:1055`),
  Q_b units J vs mW (`yelmo_thermodynamics.f90:360`), omega formula
  (`deformation.f90:383`), `t_dep_source "tracer"` (`yelmo_defs.f90:725`),
  stale OpenMP NaN warning above an active directive (`yelmo_thermodynamics.f90:705`),
  chion tsrf "as smbpal ITM" (gross refreezing is intentional;
  `surface_chion.f90:34,373`), mask long_name `z_srf > z_sl`
  (`kryos_output.f90:354`), `step_isostasy` header, `marine_shelf` set_ocean_mask
  lake code, pico T_box units, `calving_aa.f90:229,604`.

## Suggested release path

1. Pre-tag, answer-neutral except where stated: A1-A4, B5-B7, C8-C12, D13-D21,
   E33-E34, H (housekeeping). One or two commits per repo.
2. Decide now: E25 (freezing-point depth), E27 (tf_corr round trip), E28
   (K24 restart: fix or document), E29 (pc mask threshold: defer).
3. Post-release: F and G, in roughly the order F36-F38 then G43-G44.
