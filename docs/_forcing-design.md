# Forcing and climate in yelmox: a unified design

Status: proposal (2026-10-18), nothing implemented. The leading `_` keeps
quarto from rendering it. Merges and supersedes the open parts of
`_snapclim-retirement-plan.md` and `_grl-paleo-plan.md`, and takes in the LIS
work (`paper-lis-dynamics`, on 2.4-beta).

States checked: yelmox dev 0ebbd14, fesm-utils origin/dev 021f403 (the local
clone was at b8a9a3d), FesmData b420a94 (`Paleoclimate/remap.toml`,
`shared/domains.toml`), LIS `par/yelmox_LIS.nml` and `sims-log.md`
(2026-10-17).

Decided (2026-10-18): dev only, 2.4-beta frozen; indices named `idx_<name>`;
insolation on the forcing clock with an optional fixed time or orbit;
climate modifiers after the backend; a shared database with the paleoclimate
snapshots read from global products.

## 1. Use cases to cover

| Family | Domain | Climate today | Driving series | Extras |
|---|---|---|---|---|
| PD spin-up / opt (GRL, ANT) | GRL, ANT | snapclim `anom` with zeros, or snapesm `pd` | none | opt, relax |
| ISMIP7, 1pctCO2, TIPMIP, K24 | GRL, ANT | esm | esm periods | climate smb, ocean at shelf base |
| ANT paleo (spin-up, LGP) | ANT-32KM | snapclim `snap_1ind` / `fraction` | glacial index file | transient insolation, sea level |
| GRL paleo A (glacial cycles) | GRL-PAL | planned snapesm | glacial (+ millennial) index | border ice in time |
| GRL paleo B (Tabone deglac.) | GRL-PAL | planned snapesm, time-varying snapshots | none (recon in time) | Badgeley anomaly fields, border ice, LGM steady state |
| LIS (B18 M2/M3, MIS-3, S23) | LIS/NH | snapesm, 2 axes, 5 snapshots | a, b: const (steady) or NGRIP series (transient) | insolation fixed or transient, PaleoMIST init, CLIMBER-X b axis |
| North, Pyrenees, SRG | NH, PYR, SRG | snapclim | glacial index; Pyrenees orbital × 10 K | |
| Hysteresis / ramps | any | snapesm + `[tsforcing]` | tsgen ramp / feedback, kill, restart every Δf | |
| rembo GRL | GRL | rembo + snapclim ocean | `[tsforcing]` → dT_sum | |
| Bipolar | ANT + NH | snapclim both | north glacial index; OBM air anomaly | OBM ocean → shelf |

## 2. What the code does now

1. **A scalar forcing series is defined in six places**, each with its own
   reader and time convention:
   - `[tsforcing]` + `[tsgen]`: one `f_now`, three gains → `dTa/dTo/dSo`
     (`libs/kryos_forcing.f90`), driver-owned, with feedback, kill switch and
     restarts every Δf;
   - snapesm `&snap_idx_<name>`: one private tsgen per axis
     (`libs/snapesm.f90`), response variable a dummy 0;
   - snapclim `fname_at/ap/ao/as/bt/...` files;
   - OBM: `climate_air_anom` = snapclim `at` × `dTa_const` (works for
     snapclim only);
   - chion `const_insol` / `const_kabp`: a second forcing clock next to
     `ctrl.tstep_const` (the LIS insolation bug of 2026-10-17 came from this);
   - barystatic sea level from file (FastIsostasy `bsl`, on `time_rel`: fine).

   Consequence in the LIS config: `at`/`ap` and `bt`/`bp` are identical
   16-line groups, because a series is bound to one axis.
2. **The homogeneous anomaly is applied three ways**: snapclim in `anom`
   modes only, snapesm always and before the elevation transform, rembo with
   its own mapping (dT_sum, 1.3 × for annual, f_to for the ocean).
3. **snapesm gaps** (from the retirement plan, confirmed): ocean only monthly
   and only with salinity; `zs/tas/pr` required in every snapshot (so no
   ocean-only use for rembo); `dTa/dTo/dSo_const` left over; the `fraction`
   ocean rule reads the previous update's `ta_ann` (snapclim's lag
   reproduced; snapclim marks the rule "untested"); `f_stdev` inactive;
   `snapesm_end` a stub.
4. **Configuration volume**: every par file carries its own database groups
   (LIS: ~350 lines of `clim_*`), every tsgen key is required, and yelmox
   groups have no defaults file (yelmo groups do: `input/yelmo_defaults.nml`
   via `nml_read(..., defaults_file=)`, already in fesm-utils `nml`).
5. **Retirement plan, step 1 is done**: tsgen `norm_t0`/`norm_t1` (point or
   window) is on fesm-utils origin/dev and `snapesm-axes` (021f403). The plan
   found it missing because the local clone was not fetched.

The rest is in good shape and stays: kryos domain/coupler/hub, the per-step
sequence written out in the driver, `climate_out_class` as the only thing
components read, varslice with time-varying groups and online remap,
`[climate_bias]`, the esm backend (`esm_forcing.f90` untouched).

## 3. Concepts

Five nouns a user needs:

1. **Forcing clock**: `ts%time_rel`. `ctrl.tstep_method = "const"` holds it
   at `tstep_const`, `"rel"` follows the model time. Everything time-dependent
   (indices, time-varying snapshots, insolation, sea level, recon slices,
   border ice) is evaluated on it. No second clock; the only exception is an
   explicit fixed orbit (section 5).
2. **Indices**: named scalar series, `[idx_<name>]`, listed in
   `[indices] names = ...`. Owned by the driver (like `bsl`), updated once per
   step on the forcing clock with the response variable (ice volume), shared
   by every domain. Consumers refer to an index by name. An index is any
   scalar series, normalised or not (section 4).
3. **Data**: varslice groups in a database (`var_defs`, a list of files). A
   **snapshot** is a named set of field bindings, static or varying in time.
4. **Climate** = backend → modifiers → `climate_out`:
   - backend: snapshots (today snapesm; name open, section 6), `esm`,
     `rembo`;
   - modifiers, applied in listed order to `climate_out`, any backend:
     `anom` (homogeneous dTa/dTo/dSo = gain × index), `ocn_fraction`
     (ocean anomaly = f_to × mean air anomaly), `ocn_obm` (bipolar), `bias`
     (existing `[climate_bias]`); later e.g. `pr_var` (f_stdev), regional
     tapers.
5. **Prescribed ice**: `[ice_bnd]`, a data group (ice history) + a region,
   optionally scaled by an index; sets yelmo `bnd%H_ice_fix` (GRL plan).

Per step (unchanged order, one new line):

```
bsl_update(bsl, time_rel)
indices_update(idx, time_rel, var=V_ice)        ! replaces tsforcing_update
step_relax / step_optimize / step_isostasy / couple_to_yelmo
step_ice_bnd(dom, ts, idx)                       ! new, GRL-PAL only
step_icesheet / couple_yelmo_to_htopo
step_climate(dom, ts, idx)                       ! backend, then modifiers
step_surface / step_shelf
```

## 4. Indices

An index is a tsgen group with a name. What it holds depends on the method:

| Kind | tsgen | Value | Example |
|---|---|---|---|
| position on a snapshot axis | `series` with `norm_t0`/`norm_t1`, or `const` | 0 at `lo`, 1 at `hi` | glacial index, B18 a and b |
| physical series, as is | `series` without normalisation | in the file's units | a temperature anomaly [K], e.g. NGRIP or a global mean |
| generated | `const`, `ramp-*`, `sin` | any | sensitivity ramp |
| feedback | `PI42`, `H312PID`, ... | controlled by the response (ice volume) | hysteresis |
| from a coupled model (later) | not tsgen: set by the driver each step | the model's output | OBM box temperature |

A consumer applies its own gain, so a raw series in K drives `anom` with
`f_ta = 1`. An optional `units` key labels the index in the output
(`idx_<name>` in the 1D file). A series of *fields* (a 2D anomaly in time) is
data, not an index: a time-varying snapshot (section 6).

`with_kill` and `restart_every_df` stay per index; the run stops when any
index trips its kill switch. Restart: `idx_<name>_restart.nc` in the bundle,
for every index (replaces `tsgen_restart.nc` and
`snapesm_idx_<name>_restart.nc`).

## 5. Insolation and orbit

All insolation users (chion, smbpal) take an **orbit** instead of a time:

```fortran
&orbit
    method = "clock"     ! "clock": orbit at the forcing clock (Berger tables)
                         ! "time":  orbit at a fixed time [yr, model time]
                         ! "fixed": ecc, obl, per given
    time   = -21e3
    ecc    = 0.0167
    obl    = 23.44       ! [deg]
    per    = 102.7       ! [deg] longitude of perihelion
/
```

- Replaces `const_insol`/`const_kabp` in `[surface_chion]` and smbpal.
- `libs/insol/insolation.f90` already separates the two steps
  (`calc_orbital_par(time) -> ecc, obl, per`, then the daily insolation); the
  change adds an entry point taking `(ecc, obl, per)` directly.
- To share it with CLIMBER-X later, the module moves to fesm-utils
  (`insolation`, with an `orbit_class`), the orbital tables as its data
  (question 3).

## 6. Climate backends

### Snapshots (today snapesm) and esm: both are needed

| | Snapshots (snapesm) | esm (`esm_forcing.f90`) |
|---|---|---|
| Idea | climate states blended by indices | one ESM run as anomalies to an observed reference, following a protocol (ISMIP6/7) |
| Time | the forcing clock drives indices; snapshots static or slowly varying | calendar periods: reference, historical, projection; ESM reference period |
| Fields | air temperature, precipitation, ocean profiles | also smb anomaly and its elevation gradient, ocean at the shelf base with basin extrapolation, subglacial discharge |
| Extras | elevation transform, precipitation sensitivity, axes | interannual variability resampled from the ESM (`clim_var`, seeded) |
| Use | paleo, PD, sensitivity | ISMIP7, 1pctCO2, TIPMIP, K24 |

Structurally the esm atmosphere is a snapshot axis with a time-varying `hi`
(the ESM), a static `lo` (its reference period) and index 1, so the two could
converge later. Not proposed now: the smb, shelf-base ocean, variability and
discharge are protocol-specific, and the esm configs are verified against
ISMIP7.

Name: "snapesm" reads as a relative of "esm" (it was snapclim rebuilt on the
fesm-utils primitives); "snap" alone is too broad. Proposal: `snapshots`
(`[comps] climate = "snapshots"`, `libs/snapshots.f90`; groups keep the
`snap_` prefix). Question 1.

### The snapshot backend after the change

- An **axis** names two snapshots and an index:
  `&snap_axis_<name>: lo, hi, index`. A field lists its axes.
  `linear: now = ref + Σ a_k (hi_k - lo_k)`,
  `ratio: now = ref Π (1 + a_k (hi_k/lo_k - 1))` (as now, B18).
- `lo = ""`: `hi` is itself the anomaly (linear) or ratio (ratio). Covers
  Badgeley2020 (anomaly + precipitation fraction) without a zero snapshot.
- Time-varying snapshots as now (`with_time`, `extrap`): Buizert2018 as `hi`
  over a static baseline `lo`, index `const` 1.
- Fields are optional per snapshot; the backend builds only what is bound
  (ocean-only for rembo).
- Ocean: annual or monthly, salinity optional (0 when absent, as snapclim).
- Gone from `[snap]`: `dTa/dTo/dSo_const`, `f_to` (to `ocn_fraction`), the
  private tsgen indices, the `dTa/dTo/dSo` arguments of `snapesm_update`.

Config sketch, LIS B18 (today ~140 lines of `[snap*]` groups):

```fortran
&indices
    names = "a" "b"
/
&idx_a
    method = "const"
    f_min  = 0.87                 ! 0: PI, 1: LGM
/
&idx_b
    method = "series"
    series_file = "forcing/indices/kv_m3_b_c3.dat"
/
&snap
    var_defs  = "input/clim_db.nml"
    ref_name  = "pd"
    snapshots = "pd" "picontrol" "lgm" "lgmweak" "lgmstrong"
    lapse     = 0.0080 0.0065
    f_p       = 0.05
/
&snap_axis_a
    lo    = "picontrol"
    hi    = "lgm"
    index = "a"
/
&snap_axis_b
    lo    = "lgmweak"
    hi    = "lgmstrong"
    index = "b"
/
&snap_field_tas
    blend = "linear"
    axes  = "a" "b"
/
&snap_field_pr
    blend = "ratio"
    axes  = "a" "b"
/
&snap_lgm
    tas = "pmip3_lgm_mean_tas_ann" "pmip3_lgm_mean_tas_jja"
    pr  = "pmip3_lgm_mean_pr_ann"
    zs  = "pmip3_lgm_mean_zs"
/
&climate_mods
    names = "ocn_fraction"
/
&ocn_fraction
    f_to = 0.25
/
```

Unset keys come from `input/yelmox_defaults.nml`.

Pyrenees (snapclim `anom`, orbital × 10 K):

```fortran
&indices
    names = "orb"
/
&idx_orb
    method      = "series"
    series_file = "input/alpha_orbital_t2018.dat"
/
&climate_mods
    names = "anom"
/
&anom
    index = "orb"
    f_ta  = 10.0
/
```

Hysteresis ramp (today `[tsforcing]`): the same with a `ramp-*` or feedback
tsgen.

## 7. Modifiers: what changes, by example

Modifiers act on `climate_out` after the backend, in the order of
`[climate_mods] names`.

**`anom` on the snapshot backend: same result.** Today dTa is added at sea
level before the elevation transform (Γ lapse rate, β precipitation
sensitivity, z model surface); after the change, to the output:

```
today:  tas = (tsl + dTa) - Γ z     pr = prcor e^(β dTa) e^(β (tas - tsl - dTa)) = prcor e^(β dTa) e^(-β Γ z)
new:    tas = (tsl - Γ z) + dTa     pr = [prcor e^(-β Γ z)] e^(β dTa)
```

Equal up to rounding. Example: Pyrenees at 21 ka with the orbital index at
0.9: dTa = 9 K in both, the same `tas` and `pr`.

**`ocn_fraction`: no lag.** Today (snapesm, as snapclim) the ocean anomaly of
an update uses `ta_ann` of the previous update, and `now` is seeded with the
reference at load, so the first update after a cold start *and after every
restart* sees no anomaly. LIS LGM steady state (a = 0.87, b = 0.48,
`dt_clim` 10 yr, f_to 0.25): today the ocean stays at the present-day
reference for the first 10 yr after the cold start and after each restart,
then follows the air anomaly one update late; after the change it follows
from the first update. In a transient the lag is one `dt_clim` (10 yr): no
visible difference. Order: with `names = "anom" "ocn_fraction"` the ocean
also sees the homogeneous anomaly (as today, one update late).

**What the fraction rule averages (question 4).** Today: the mean of
`ta_ann - ta_ann_ref` over *all* cells of the climate grid, with `ta_ann` at
the current model surface and `ta_ann_ref` at the reference snapshot's
surface. So it includes:

- the lapse-rate cooling of a growing ice sheet: for a LIS of ~2 km over
  ~30% of the climate grid, roughly 0.3 × 2000 m × 8 K/km ≈ 5 K of mean
  cooling, ≈ 1.2 K more ocean cooling at f_to 0.25 (estimate);
- at present day, any difference between the model surface and the
  reference orography (e.g. ERA5 at its own resolution): a non-zero ocean
  anomaly with no forcing.

Proposal: average the sea-level anomaly `tsl_ann - tsl_ann_ref` (the
forcing only), optionally over a region (`ocn_fraction.region`, e.g. ocean
cells). This changes the results of every `fraction` config.

## 8. Bipolar: the ocean box model

Today, per hemisphere, the OBM box temperatures are set from the domain's
air anomaly dTa (snapclim `at` × `dTa_const`) with the regional
amplification factors (Zickfeld et al., 2004): `thetan = thetan_init + dTa`
(north), `thetat = thetat_init + dTa ampl_tropics/ampl_north`, vapor fluxes
likewise; the south does the same with `ampl_south` and overwrites `thetat`
when both are active. Switches `atm2obm`, `ism2obm`, `obm2ism` in `[ctrl]`.
`obm2ism` overwrites the domain's `to_ann` (all depths) with the box
temperature, outside the climate step.

Proposal, with a global anomaly as the single input:

```fortran
&obm_forcing
    index = "gmt"          ! global mean air anomaly [K] (an index)
    gain  = 1.0
/
```

- `atm2obm`: every box from the global anomaly, `theta_i = theta_i_init +
  ampl_i × dT_glob` (north, tropics, south), vapor fluxes as now; the
  tropical box is no longer set twice.
- The ice-sheet domains take the same index with their own regional gain
  (`anom` with `f_ta = ampl_north`, or a snapshot axis): explicit in the par
  file.
- `obm2ism` becomes the modifier `ocn_obm` (mode `absolute`, as now, or
  `anomaly`: reference ocean + box anomaly, keeping the pattern); the box
  temperatures enter as model-set indices (`obm_tn`, `obm_ts`, section 4).
- Switches stay; the order of the exchange is unchanged.

## 9. Shared database

`input/clim_db.nml` (one file, in yelmox): one group per product field,
named after the product. Configs list databases: `var_defs =
"input/clim_db.nml" "<par file>"`, the par file holding only run-specific
groups (e.g. the LIS NGRIP experiments in the paper repo). The esm backend
keeps its own par_file for now.

- **Present-day references stay regional** (high resolution), with
  `{domain}/{grid_name}` substitution: MAR v3.14, ERA5, RACMO, ORAS5, ISMIP6
  ocean.
- **Paleoclimate snapshots from the global products**, remapped online
  (bilinear, map cached in `maps/`), one group set for every domain:

| Product (FesmData) | Global grids | Groups |
|---|---|---|
| `Paleoclimate-PMIP3-<exp>-mean` (lgm, midHolocene, piControl) | GLOBAL-0.5/1/2/5DEG | `pmip3_<exp>_mean_{tas_ann,tas_jja,tas_djf,pr_ann,zs}` |
| `Paleoclimate-PMIP4-<exp>-mean` | GLOBAL-0.5/1/2/5DEG | `pmip4_<exp>_mean_{tas_mon,pr_mon,zs}` (or ann/jja/djf) |
| `Paleoclimate-CLIMBER3a-<sim>` (present, lgm_1p7weak, lgm_1p7strong), `-ocean` | atmosphere 0.5/1/2/5DEG, ocean 1/2/5DEG | `c3a_<sim>_{tas_ann,tas_jja,pr_ann,zs}`, `c3a_<sim>_to` |
| `Paleoclimate-CLIMBERX-Willeit2025-<state>` (piControl, interstadial, stadial, heinrich), `-ocean` | GLOBAL-1/2/5DEG (source 5°) | `cx_<state>_{tas_mon,pr_mon,zs}`, `cx_<state>_{to,so}` |

- Global grid: GLOBAL-1DEG for all (the finest common to the four
  families; the sources are about 1° or coarser, CLIMBER-X 5°), so every
  snapshot of a run is remapped with the same kernel. Whether the regional
  `clim_sigma` smoothing of PMIP3 (50 km) is still needed is checked in the
  climate-only comparison (section 13); `clim_sigma` stays an option.
- Southern summer: the Antarctic snapshots need DJF (`t2m_djf`, PMIP3/4) or
  monthly fields (PMIP4, CLIMBER-X); CLIMBER-3α has JJA only.
- The single PMIP models (`-<model>-<exp>`) and `-stdev` are added in the
  same form when needed.

## 10. Use cases in the new form

| Case | Backend | Indices | Modifiers | Other |
|---|---|---|---|---|
| PD (snapclim `const`/`anom` zeros) | snapshots `pd` only | none | `ocn_fraction` | |
| ISMIP7 / 1pctCO2 / TIPMIP / K24 | esm | none | `bias` (opt.) | unchanged |
| ANT paleo | snapshots pd/picontrol/lgm, axis a | a = series glacial index | `ocn_fraction` | orbit `clock` |
| North / LIS legacy `snap_1ind` | as ANT paleo | a | `ocn_fraction` | |
| LIS B18 | snapshots, axes a, b (b: CLIMBER-3α or CLIMBER-X) | a, b: const or NGRIP series | `ocn_fraction` | orbit `clock` or `time` |
| Pyrenees | snapshots `pd` | orb | `anom` f_ta 10 | |
| SRG (`const`) | snapshots `pd` | none | none | |
| GRL paleo A | snapshots, axis a (+ b) | glacial (+ millennial) | `ocn_fraction` or ocean axis | `[ice_bnd]` |
| GRL paleo B | snapshots, Buizert `hi` over baseline `lo`; Badgeley `lo=""` | const 1 | ocean: question 7 | `[ice_bnd]`, LGM steady = `tstep_method const` at -21 ka |
| rembo GRL | rembo + snapshots ocean-only | ramp index | rembo reads the `anom` gains itself | |
| Bipolar | snapshots per domain | gmt (+ glacial) | `anom` per domain, `ocn_obm` | `[obm_forcing]` |
| Ramps / hysteresis | snapshots | feedback index | `anom` | kill, restart every Δf |

## 11. Extension points

To test a new idea, one of:

| New... | Where | Touches |
|---|---|---|
| time series shape (window, taper, sum of series) | tsgen (fesm-utils) | nothing in yelmox |
| spatial pattern / snapshot | `input/clim_db.nml` group + snapshot or axis | nothing in code |
| operation on the climate (variability, regional taper, ...) | one subroutine + one `case` in `libs/climate_mods.f90` | nothing else |
| orbit | `[orbit]` | nothing in code |
| climate model | backend in `yelmox_climate.f90` | fills `climate_out` |
| prescribed ice / boundary | `[ice_bnd]` data + region | |

Region-specific behaviour stays opt-in by key, never `select case(domain)`.

## 12. Code changes

| File | Change |
|---|---|
| `libs/kryos_forcing.f90` | `tsforcing_class` → `indices_class`: named tsgen array (+ model-set entries), `indices_init/update/value(name)/restart_*`, kill and Δf restarts per index; `idx_<name>_restart.nc` |
| `libs/snapesm.f90` (→ `snapshots.f90`) | axes reference indices by name; `lo=""`; optional fields; annual/no-salinity ocean; drop dTa args, `*_const`, `f_to`, fraction lag; `snapesm_end` |
| `libs/climate_mods.f90` (new) | `anom`, `ocn_fraction`, `ocn_obm`, `bias` (moved from `yelmox_climate.f90`) |
| `libs/yelmox_climate.f90` | backends snapshots/esm/rembo; `climate_update(..., idx)` then `climate_mods_apply`; snapclim cases removed; rembo ocean via snapshots |
| `yelmox_bipolar/obm_coupling.f90` | `[obm_forcing]` global anomaly → boxes; `obm2ism` → `ocn_obm`; `climate_air_anom` removed |
| `libs/insol/` → fesm-utils `insolation` | `orbit_class`, insolation from `(ecc, obl, per)`; chion/smbpal take `[orbit]`, `const_insol/const_kabp` removed |
| `libs/ice_bnd.f90` (new) + yelmo `bnd%H_ice_fix` | GRL paleo border ice |
| `input/clim_db.nml`, `input/yelmox_defaults.nml` (new) | shared database (section 9); defaults for yelmox groups, read with `defaults_file` |
| `scripts/check_par_nml.py` | also check yelmox groups against `yelmox_defaults.nml` |
| removed | `libs/snapclim.f90`, its tests, `libs/obm/atm2ocean.f90`, `input/greenland_snp.nml` |

## 13. Retiring snapclim with confidence

1. Coverage: every row of the retirement plan's usage table has a converted
   config (section 10); `check_par_nml.py` clean.
2. Climate-only comparison, before deleting anything: each converted config
   and its snapclim original with `with_ice_sheet = False`, `tstep_method =
   "rel"`, `dtt = dt_clim = 1000` over the full forcing range (-130 ka to 0
   for paleo, 0 for PD): `snap.nc` of both, maps and domain means of
   `t2m_ann`, `pr_ann`, `to_ann` per time. Minutes per config; shows where
   and why they differ (dropped `f_corr`, fraction rule, global vs regional
   remap of the snapshots).
3. 10-yr coupled runs on Levante per family (as planned).
4. Delete snapclim in its own commit, after 2-3 pass.

## 14. Phases and commits

All on dev; 2.4-beta stays frozen.

| # | Repo | Content | Check |
|---|---|---|---|
| 0 | yelmox | this doc; retirement plan step 1 marked done | |
| 1 | yelmox | indices registry; snapesm axes on indices; `[tsforcing]` → `anom` modifier; per-index restarts | GRL snapesm, LIS steady: identical up to rounding; ramp restart continuity |
| 2 | yelmox | `climate_mods`: `ocn_fraction` (no lag, question 4), `bias` moved | expected diffs only |
| 3 | yelmox | snapshot backend: optional fields, annual/no-salt ocean, `lo=""`; rembo ocean on it | |
| 4 | fesm-utils + yelmox | insolation with `orbit_class` in fesm-utils; `[orbit]` in chion/smbpal | `clock` identical to transient insolation today, `time` to `const_insol` |
| 5 | yelmox | `clim_db.nml` (global paleoclimate), `yelmox_defaults.nml`, checker | configs shorter |
| 6 | yelmox | config conversions per family (retirement plan order); bipolar with `[obm_forcing]`, `ocn_obm` | section 13 steps 1-3 |
| 7 | yelmox | remove snapclim; rename snapesm (question 1); docs, CHANGELOG | |
| 8 | fesm-utils | tsgen windowed modifiers (replaces `f_hol`), RNG state in restart | |
| 9 | yelmo + yelmox | `H_ice_fix`, `[ice_bnd]`, GRL-PAL base config, GRL paleo A then B | |

LIS moves to the new keys after phases 1-5 (paper repo par files, on dev).

## 15. Decisions needed

1. Backend name: keep `snapesm`, or `snapshots` (proposed)?
2. Indices set by a coupled model (OBM box temperatures) in the same
   registry, or kept separate?
3. Insolation module to fesm-utils (shared with CLIMBER-X later), or stays in
   yelmox for now?
4. `ocn_fraction`: average the sea-level anomaly (forcing only, proposed) or
   keep today's surface-temperature mean over the whole grid (includes the
   ice sheet's own elevation cooling)? Optional region?
5. OBM: one global anomaly into all boxes via the amplification factors
   (proposed), or keep the per-hemisphere input? `ocn_obm` absolute (as now)
   or anomaly?
6. Global grid of the paleoclimate snapshots: GLOBAL-1DEG for all
   (proposed)? Antarctic summer from DJF/monthly fields where available?
7. GRL paleo B ocean: `ocn_fraction`, CLIMBER snapshots, or a TraCE product
   (GRL plan Q9; its other questions stay open).
