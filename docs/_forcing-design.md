# Forcing and climate in yelmox: a unified design

Status: proposal (2026-10-18, revised after a code review the same day),
nothing implemented. The leading `_` keeps quarto from rendering it. Merges
and supersedes the open parts of `_snapclim-retirement-plan.md` and
`_grl-paleo-plan.md`, and takes in the LIS work (`paper-lis-dynamics`, on
2.4-beta).

Review outcome (2026-10-18): the architecture stands; the core decisions match
seams already in the code (`climate_out_class` as the only consumer, the
`[climate_bias]` apply step as the modifier pattern, varslice `subs` for path
keywords, `nml` `defaults_file`/`defaults_group`). Revisions: modifiers in a
fixed order with `active` keys, not an ordered list (section 7); a written
contract per modifier (section 7); the OBM box temperatures passed directly,
not registry indices (sections 4, 8); a tsgen defaults passthrough in
fesm-utils (section 14); phase 1 split into a behaviour-preserving and a
behaviour-changing commit; the climate-only comparison as a permanent script
(section 13); esm on calendar time named as the clock exception (section 3).

States checked: yelmox dev 0ebbd14, fesm-utils origin/dev 021f403 (the local
clone was at b8a9a3d), FesmData b420a94 (`Paleoclimate/remap.toml`,
`shared/domains.toml`), LIS `par/yelmox_LIS.nml` and `sims-log.md`
(2026-10-17), fesmc/insol main 27e767e, configme `packages/insol.toml`.

Decided (2026-10-18): dev only, 2.4-beta frozen; backend renamed
`snapshots`; indices named `idx_<name>`, any scalar series; insolation from
the fesmc/insol library (a configme dependency), on the forcing clock with an
optional fixed time or orbit; climate modifiers after the backend, in a
fixed order, each an optional group with `active`; `ocn_fraction` on the
sea-level anomaly; OBM driven by one global anomaly, its box temperatures
passed to `ocn_obm` by the bipolar bridge (no model-set indices until a
second user exists); a shared database with path keywords, the paleoclimate
snapshots read from global products (GLOBAL-1DEG).

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
   16-line groups, because a series is bound to one axis. (This is a
   snapclim limitation; in snapesm two fields can already share one index.
   The axis layer of section 6 is justified by driver ownership of the index,
   needed by modifiers, `[ice_bnd]` and multi-domain drivers, not by the
   duplication.)
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
   via `nml_read(..., defaults_file=)`, already in fesm-utils `nml`, with a
   `defaults_group` argument that maps a named group such as
   `snap_field_tas` onto a template group `snap_field`). `tsgen_init` reads
   its keys with plain `nml_read`, so `[idx_<name>]` groups stay fully
   required until tsgen passes `defaults_file`/`defaults_group` through
   (phase 8).
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
   border ice) is evaluated on it. No second clock. Two explicit exceptions:
   a fixed orbit (section 5), and the esm backend, which follows calendar
   time (`ts%time` against `calendar_ref`) because its periods are calendar
   years; esm configs run with `tstep_method = "rel"` and `calendar = True`,
   so the two clocks coincide.
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
   - modifiers, applied to `climate_out` in a fixed order, each an optional
     group with `active` (as `[climate_bias]` today): `anom` (homogeneous
     dTa/dTo/dSo = gain × index), `ocn_fraction` (ocean anomaly = f_to ×
     mean sea-level air anomaly), `ocn_obm` (bipolar), `taper` (regional
     weight towards the reference), `bias` (existing `[climate_bias]`);
     later e.g. `pr_var` (f_stdev). Each modifier states which backends it
     accepts (section 7); modifier groups take the domain suffix
     (`anom<sfx>`), as `climate_bias<sfx>` does.
5. **Prescribed ice**: `[ice_bnd]`, a data group (ice history) + a region,
   optionally scaled by an index; writes yelmo `bnd%H_ice_ref` in the
   fixed/relaxation strip (`border_width`, GRL plan).

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

Every index is a tsgen series. A scalar set by another model each step (OBM
box temperature, later a coupled climate model's global mean) is coupler
state, not an index: it is passed to its consumer by the driver that owns
the model (section 8). Model-set entries in the registry would need a second
kind, their own restart path and output for one consumer today; they are
added when a second user exists.

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
- **Library: fesmc/insol, branch `main`.** yelmox still vendors an older copy under
  `libs/insol/` (module state, the same Laskar tables);
  fesmc/insol is its successor: a standalone library (no netCDF, no
  fesm-utils), state object `insol_class` / `orbit_class`, Laskar et al.
  (2004) elements, daily, diurnal and daylength routines, already a configme
  package (`configme install insol`, `libinsol.a`), but listed as
  "standalone" because yelmox vendors the source.
- Change: insol becomes a yelmox dependency (`.configme/manifest.toml`
  `deps += "insol:<branch>"`, the configme package entry no longer
  "standalone"; `config/Makefile_yelmox.mk` links `libinsol.a`, includes its
  modules and drops the `interp1D`/`insolation` rules); `libs/insol/` is
  removed; chion's host insolation (`surface_chion.f90`) and smbpal call
  insol with an `insol_class` initialised once.
- New functionality goes into insol: `orbit_class` from fixed `(ecc, obl,
  per)` and an insolation call taking an orbit instead of a time (all calls
  take `time_bp` today). CLIMBER-X can use the same library later.
- Check: insol has changed since the copy was vendored (precision and solar
  constant unified, selectable mean-cosine convention): compare the daily
  insolation of both at 0, 21, 40 and 125 ka before switching. A good
  match is enough: small differences (interpolation etc.) are accepted.

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
fesm-utils primitives); "snap" alone is too broad. Decided: `snapshots`
(`[comps] climate = "snapshots"`, `libs/snapshots.f90`; groups keep the
`snap_` prefix).

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
&ocn_fraction
    active = True
    f_to   = 0.25
/
```

Unset keys come from `input/yelmox_defaults.nml`; named groups
(`snap_field_<name>`, `snap_<name>`, `snap_axis_<name>`, `idx_<name>`) fall
back to template groups there (`snap_field`, `snap_snapshot`, `snap_axis`,
`idx`) via `defaults_group`, so `blend = "linear"` and `axes = ""` need not
be written for a field that only follows the reference.

Pyrenees (snapclim `anom`, orbital × 10 K):

```fortran
&indices
    names = "orb"
/
&idx_orb
    method      = "series"
    series_file = "input/alpha_orbital_t2018.dat"
/
&anom
    active = True
    index  = "orb"
    f_ta   = 10.0
/
```

Hysteresis ramp (today `[tsforcing]`): the same with a `ramp-*` or feedback
tsgen.

## 7. Modifiers: order, contract, and what changes

Modifiers act on `climate_out` after the backend, in a **fixed order**:
`anom`, `ocn_fraction`, `ocn_obm`, `taper`, `bias`. Each is an optional
namelist group with `active` (default off), read like `[climate_bias]`
today. An ordered `names` list was considered and dropped: only one order is
physically meaningful, a fixed order removes a mis-ordering failure mode and
one group, and the par checker can validate each group against the defaults
file. `climate_mods_apply` is one subroutine that calls each active modifier
in that order; a new modifier is one subroutine and one line there.

### Contract per modifier

What each one touches in `climate_out%now`, what it needs, and when it
refuses (an error at init, never a silent partial application):

| Modifier | Touches | Needs | Refuses |
|---|---|---|---|
| `anom` | `tas`, `ta_ann`, `tsl_ann` (+dTa); `pr`, `pr_ann` (× e^(β dTa)); ocean: `to_ann`/`so_ann` (+dTo/+dSo) when profiles, else `T_shlf`/`S_shlf` and `dT_shlf`/`dS_shlf` | an index; `beta_p` | a backend that supplies the smb (`has_smb`: esm with `use_smb`, rembo), since a post-hoc air anomaly cannot reach that smb; rembo reads `[anom]` itself instead (below) |
| `ocn_fraction` | `to_ann` = `ref%to_ann` + f_to × mean(`tsl_ann` − `ref%tsl_ann`) over `region` | profile ocean and `tsl_ann` (snapshot backend) | `has_ocn_shelf` (esm); a backend without `tsl_ann` |
| `ocn_obm` | `to_ann` (all depths): `absolute` = box value, `anomaly` = `ref%to_ann` + (box − box_init) | profile ocean; the box value, passed by the bipolar bridge (section 8) | `has_ocn_shelf` |
| `taper` | atmosphere: `now = ref + w (now − ref)` for `tas`, `pr` and their annual means | a weight field `w` from a var_defs group (static) | nothing |
| `bias` | as today: `tas`, `ta_ann`, `tsl_ann` (+dtas), `pr`, `pr_ann` (× fpr) | var_defs groups | rembo (annual fields only) |

Everything downstream is unchanged: the elevation dependence (`lapse`,
`beta_p`, `z_srf`) that moves the atmosphere to the surface grid is additive
and multiplicative in the same way, and `smb_simple` reads the modified
`tsl_ann`.

**rembo is the one backend where the homogeneous anomaly is an input, not a
modifier**: REMBO needs `dT_sum` before it runs. The rembo backend reads the
same `[anom]` group (index and `f_ta` → `dT_sum`, `1.3 ×` annual, `f_to` →
ocean) and `climate_mods_apply` skips `anom` for it. Documented as the
exception; keeping the anomaly a backend input everywhere (today's design)
was rejected because it keeps dTa arguments in every backend and blocks the
esm case.

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
visible difference. Order: `anom` runs before `ocn_fraction`, so the ocean
also sees the homogeneous anomaly (as today, one update late).

**What the fraction rule averages.** Today: the mean of
`ta_ann - ta_ann_ref` over *all* cells of the climate grid, with `ta_ann` at
the current model surface and `ta_ann_ref` at the reference snapshot's
surface. So it includes:

- the lapse-rate cooling of a growing ice sheet: for a LIS of ~2 km over
  ~30% of the climate grid, roughly 0.3 × 2000 m × 8 K/km ≈ 5 K of mean
  cooling, ≈ 1.2 K more ocean cooling at f_to 0.25 (estimate);
- at present day, any difference between the model surface and the
  reference orography (e.g. ERA5 at its own resolution): a non-zero ocean
  anomaly with no forcing.

Decided: average the sea-level anomaly `tsl_ann - tsl_ann_ref` (the
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
- `obm2ism` becomes the modifier `ocn_obm` (mode `absolute` by default, as now, or
  `anomaly`: reference ocean + box anomaly, keeping the pattern). The box
  temperature is **passed directly** by the bipolar bridge:
  `coupling_obm2ism` calls `climate_mod_ocn_obm(dom%clim, theta, mode)` in
  `obm_exchange`, after both domains' climate steps and before their shelf
  steps, exactly where `obm2ism` runs today. So `ocn_obm` is the one
  modifier applied outside `climate_mods_apply`, by the driver that owns the
  box model; it keeps the same contract (section 7). No registry entry, no
  extra restart path: the box state is in `obm_restart.nc` already.
- Switches stay; the order of the exchange is unchanged.

## 9. Shared database

`input/clim_db.nml` (one file, in yelmox): one group per product field.
Configs list databases: `var_defs = "input/clim_db.nml"`; the par file first,
then the var_defs files in order (phase 6), the par file holding only run-specific groups (e.g. the LIS NGRIP experiments in
the paper repo). The esm backend keeps its own par_file for now.

**No group per grid or domain.** Paths use varslice keywords (built in:
`{domain}`, `{grid_name}`; extra `{key}` substitutions are already supported,
`varslice_init_nml(..., subs=)`, also in `grid_src`). Two extra keys:

- `{grid_global}`: the global grid of the paleoclimate products, set once in
  `[snap] grid_global` (default `GLOBAL-1DEG`);
- `{exp}`: the experiment/state of a snapshot, set in its group
  (`&snap_lgm exp = "lgm"`), so one group serves every experiment of a
  product family.

```fortran
&era5_1981_2010_tas                  ! regional reference, any domain and grid
    filename = "fesmdata/{domain}/{grid_name}/{grid_name}_Atmosphere-ERA5-1981-2010.nc"
    name     = "t2m"
    time     = 0 0 0 0 12
/
&pmip3_mean_tas_ann                  ! any PMIP3 experiment, any domain
    filename = "fesmdata/Global/{grid_global}/{grid_global}_Paleoclimate-PMIP3-{exp}-mean.nc"
    name     = "t2m_ann"
    remap    = "bilin"
    grid_src = "{grid_global}"
/
&snap_lgm
    exp = "lgm"
    tas = "pmip3_mean_tas_ann" "pmip3_mean_tas_jja"
    pr  = "pmip3_mean_pr_ann"
    zs  = "pmip3_mean_zs"
/
```

- **Present-day references stay regional** (high resolution), with
  `{domain}/{grid_name}`: MAR v3.14, ERA5, RACMO, ORAS5, ISMIP6 ocean.
- **Paleoclimate snapshots from the global products**, remapped online
  (bilinear, map cached in `maps/`):

| Product (FesmData) | Global grids | `{exp}` | Groups |
|---|---|---|---|
| `Paleoclimate-PMIP3-{exp}-mean` | GLOBAL-0.5/1/2/5DEG | lgm, midHolocene, piControl | `pmip3_mean_{tas_ann,tas_jja,tas_djf,pr_ann,zs}` |
| `Paleoclimate-PMIP4-{exp}-mean` | GLOBAL-0.5/1/2/5DEG | lgm, midHolocene, piControl | `pmip4_mean_{tas_mon,pr_mon,zs}` (and ann/jja/djf) |
| `Paleoclimate-CLIMBER3a-{exp}`, `-{exp}-ocean` | atmosphere 0.5/1/2/5DEG, ocean 1/2/5DEG | present, lgm_1p7weak, lgm_1p7strong | `c3a_{tas_ann,tas_jja,pr_ann,zs}`, `c3a_to` |
| `Paleoclimate-CLIMBERX-Willeit2025-{exp}`, `-{exp}-ocean` | GLOBAL-1/2/5DEG (source 5°) | piControl, interstadial, stadial, heinrich | `cx_{tas_mon,pr_mon,zs}`, `cx_{to,so}` |

- GLOBAL-1DEG for all (the finest common to the four families; the sources
  are about 1° or coarser, CLIMBER-3α 22.5° × 7.5°, CLIMBER-X 5°), so every
  snapshot of a run is remapped with the same kernel. Whether the regional
  `clim_sigma` smoothing of PMIP3 (50 km) is still needed is checked in the
  climate-only comparison (section 13); `clim_sigma` stays an option.
- The single PMIP models (`-<model>-<exp>`, with `{model}`) and `-stdev`
  are added in the same form when needed.

**Seasons.** CLIMBER-3α has annual means (`t2m_ann`, `pr_ann`) and JJA
(`t2m_jja`) only: no DJF, no monthly fields, neither in FesmData nor in the
original files (`fesmc/Montoya2008-CLIMBER3a`, `TS_ANN`, `TS_JJA`, `PRC_ANN`).
PMIP3 has ann, JJA, DJF; PMIP4 and CLIMBER-X are monthly. Today the
two-field synthesis (`[ann, sum]` → 12 months) takes the second field as the
hemisphere's summer. Decided (tentatively): a snapshot names the season of its second
field (`tas_season = "jja"` or `"djf"`, default: the summer of the
hemisphere), and the synthesis builds the cycle from the annual mean and that
season. CLIMBER-3α then works in the south too (JJA as the winter),
PMIP3 takes DJF in the south as now.

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
| GRL paleo B | snapshots, Buizert `hi` over baseline `lo`; Badgeley `lo=""` | const 1 | ocean: any of `ocn_fraction`, CLIMBER ocean snapshots, TraCE-21ka (selectable) | `[ice_bnd]`, LGM steady = `tstep_method const` at -21 ka |
| rembo GRL | rembo + snapshots ocean-only | ramp index | rembo reads `[anom]` itself (section 7) | |
| Bipolar | snapshots per domain | gmt (+ glacial) | `anom<sfx>` per domain, `ocn_obm` via the bridge | `[obm_forcing]` |
| Ramps / hysteresis | snapshots | feedback index | `anom` | kill, restart every Δf |

## 11. Extension points

To test a new idea, one of:

| New... | Where | Touches |
|---|---|---|
| time series shape (window, taper, sum of series) | tsgen (fesm-utils) | nothing in yelmox |
| spatial pattern / snapshot | `input/clim_db.nml` group + snapshot or axis | nothing in code |
| operation on the climate (variability, ...) | one subroutine + one call in `climate_mods_apply` (`libs/climate_mods.f90`), its group in `yelmox_defaults.nml`, a contract row in section 7 | nothing else |
| orbit | `[orbit]` | nothing in code |
| climate model | backend in `yelmox_climate.f90` | fills `climate_out` |
| prescribed ice / boundary | `[ice_bnd]` data + region | |

Region-specific behaviour stays opt-in by key, never `select case(domain)`.

## 12. Code changes

| File | Change |
|---|---|
| `libs/kryos_forcing.f90` | `tsforcing_class` → `indices_class`: named tsgen array, `indices_init/update/value(name)/restart_*`, kill and Δf restarts per index; `idx_<name>_restart.nc` |
| `libs/snapesm.f90` (→ `snapshots.f90`) | axes reference indices by name; `lo=""`; optional fields; annual/no-salinity ocean; drop dTa args, `*_const`, `f_to`, fraction lag; `snapesm_end` |
| `libs/climate_mods.f90` (new) | `anom`, `ocn_fraction`, `ocn_obm`, `taper`, `bias` (moved from `yelmox_climate.f90`); `climate_mods_init` (reads the groups, checks each contract against the backend), `climate_mods_apply` (fixed order), `climate_mod_ocn_obm` (public, for the bridge) |
| `libs/yelmox_climate.f90` | backends snapshots/esm/rembo; `climate_update(..., idx)` then `climate_mods_apply`; snapclim cases removed; rembo ocean via snapshots; rembo reads `[anom]` |
| `yelmox_bipolar/obm_coupling.f90` | `[obm_forcing]` global anomaly → boxes; `coupling_obm2ism` → `climate_mod_ocn_obm`; `climate_air_anom` removed |
| `fesm-utils` `tsgen.f90` | `tsgen_init(..., defaults_file, defaults_group)` passthrough, so `[idx_<name>]` keys can default |
| `libs/insol/` → fesmc/insol (configme dependency) | fixed `orbit_class` and orbit-based calls added in insol; chion/smbpal take `[orbit]`, `const_insol/const_kabp` removed; `libs/insol/` and its Makefile rules removed |
| `libs/ice_bnd.f90` (new) + yelmo masks `border_width` | GRL paleo border ice into `H_ice_ref` |
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
   remap of the snapshots). **Permanent**, not a one-off:
   `scripts/climate-compare/` (a runner that sets the overrides and runs two
   par files, a Julia script that plots the maps and means and prints the
   per-time differences). It becomes the regression check for every later
   modifier and backend change, and it is where the open data questions are
   answered: whether `clim_sigma` is still needed after the global remap
   (section 9), and the I/O cost per update of the time-varying GRL
   snapshots (two bracketing slices of a 1.5 GB file; from the GRL plan).
3. 10-yr coupled runs on Levante per family (as planned).
4. Delete snapclim in its own commit, after 2-3 pass.

## 14. Phases and commits

All on dev; 2.4-beta stays frozen.

| # | Repo | Content | Check |
|---|---|---|---|
| 0 | yelmox | this doc; retirement plan step 1 marked done; `scripts/climate-compare/` | the script runs on today's snapclim vs snapesm GRL pair |
| 1a | yelmox | indices registry; snapesm axes on registry indices; per-index restarts; `[tsforcing]` still feeds the backends' dTa args | behaviour-preserving: GRL snapesm, LIS steady, ramp restart continuity bit-for-bit (or up to rounding) |
| 1b | yelmox | `climate_mods` with `anom` (fixed order, `active`); `[tsforcing]` and the dTa args removed; rembo reads `[anom]`; contract checks at init | the anomaly crosses the elevation transform: equal up to rounding (section 7); esm + `anom` + `use_smb` refused |
| 2 | yelmox | `ocn_fraction` (no lag, sea-level anomaly, `region`), `bias` moved | expected diffs only |
| 3 | yelmox | snapshot backend: optional fields, annual/no-salt ocean, `lo=""`; rembo ocean on it | |
| 4 | insol + configme + yelmox | insol: fixed orbit; configme: insol a yelmox dependency; yelmox: link libinsol, `[orbit]` in chion/smbpal, `libs/insol/` removed | insolation tables vs the vendored copy (section 5); `clock` = transient insolation today, `time` = `const_insol` |
| 5 | fesm-utils + yelmox | fesm-utils: `tsgen_init` defaults passthrough; yelmox: `clim_db.nml` (global paleoclimate), `yelmox_defaults.nml` with template groups, checker | configs shorter; every `idx_<name>` group can omit unused keys |
| 6 | yelmox | config conversions per family (retirement plan order); bipolar with `[obm_forcing]`, `ocn_obm` via the bridge | section 13 steps 1-3 |
| 7 | yelmox | remove snapclim; rename snapesm → snapshots; docs, CHANGELOG | |
| 8 | fesm-utils | tsgen windowed modifiers (replaces `f_hol`), RNG state in restart | |
| 9 | yelmo + yelmox | `border_width`, `[ice_bnd]`, `taper` modifier, GRL-PAL base config, GRL paleo A then B | |
| 9b | FesmData | TraCE-21ka ocean product (time-varying `to`, `so`; not in v2, v1 had `OCN-RECON-TRACE21ka`), global, for the GRL paleo B ocean | |

LIS moves to the new keys after phases 1-5 (paper repo par files, on dev).

## 15. Decisions needed

All answered (2026-10-18): name `snapshots`; indices are tsgen series only,
the OBM box temperature is passed by the bridge (revised from "model-set
indices in the same registry" after review); modifiers in a fixed order with
`active` keys (revised from an ordered `names` list); insolation from
fesmc/insol `main` as a dependency; `ocn_fraction` on the sea-level anomaly;
OBM from one global anomaly; GLOBAL-1DEG with path keywords; `tas_season`
for two-field snapshots (tentative); `ocn_obm` absolute by default,
`anomaly` as an option; GRL paleo B ocean: all three options
(`ocn_fraction`, CLIMBER ocean snapshots, a TraCE-21ka product) selectable,
for comparison.

insol need not reproduce the vendored copy exactly (a good match is
enough). GRL paleo decisions: `_grl-paleo-plan.md`.
