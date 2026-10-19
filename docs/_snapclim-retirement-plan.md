# Retiring snapclim: snapesm + tsgen for all snapshot climates

Status: design draft (2026-10-17), not started. The leading `_` keeps quarto
from rendering it. Superseded where they differ by `_forcing-design.md`
(2026-10-18).

States checked: yelmox dev ed2a6b5, fesm-utils dev b8a9a3d, FesmData
products on Levante (`/work/ba1442/fesmdata`, paleoclimate-v2.6.0).

## Goal

- Remove `libs/snapclim.f90`. Every snapshot climate runs on `snapesm`
  (`[comps] climate = "snapesm"`), and every spatially homogeneous anomaly
  comes from the transient forcing (`[tsforcing]` + `[tsgen]`).
- The Montoya and Levermann (2008) CLIMBER-3α snapshots are read from the
  FesmData global products and remapped at load, so one database group works
  for any domain.
- No exact match with snapclim is required: 10-yr Levante runs of each
  converted config are the check.

Decided so far: anomalies via tsgen (not a snapesm index); bipolar is
converted in the same change.

## Current use of snapclim

| Namelist | atm_type / ocn_type | Notes |
|---|---|---|
| `yelmox_Antarctica`, `yelmox_pd_Antarctica`, `yelmox_pd_Greenland` | anom / anom | `fname_at` = zeros, `dTa_const` = 0: in effect the reference |
| `yelmox_Pyrenees` | anom / anom | `alpha_orbital_t2018.dat` × `dTa_const` = 10 K |
| `yelmox_Antarctica_paleo_spinup` | anom / fraction | |
| `yelmox_Antarctica_paleo_lgp`, `yelmox_LIS`, `yelmox_North` | snap_1ind / fraction | PMIP3 piControl -> lgm |
| `yelmox_Greenland`, `yelmox_Greenland_snapclim` | snap_1ind_new / fraction | snapesm port: `yelmox_Greenland_snapesm.nml` |
| `yelmox_SRG` | const / const | |
| `yelmox_rembo_Greenland` | snap_1ind_new / const | rembo atmosphere, snapclim ocean only |
| `yelmox_bipolar_Bipolar` | south anom / anom (zeros); north snap_1ind_new / fraction | OBM reads `climate_air_anom` |

Unused modes, dropped: `hybrid`, `recon`, `snap_2ind*`, `*_abs`,
`snap_1ind_miocene`.

With the current types no config reads a Montoya2008 file: `clim3` is loaded
only by `snap_2ind*`/`miocene`, the snapshot oceans only by `snap_*` ocean
types, and `clim1`/`clim2` already point to PMIP3.

## Mapping onto snapesm

| snapclim | snapesm |
|---|---|
| `const` | the `pd` snapshot only, no index |
| `anom` | the `pd` snapshot only; anomaly from `[tsforcing]` (snapesm adds dTa/dTo/dSo on top in every mode) |
| `snap_1ind`, `snap_1ind_new` | snapshots `pd`, `picontrol`, `lgm`; indices `at`, `ap` (picontrol -> lgm), normalised in tsgen |
| `fraction` (ocean) | field `to` with `blend = "fraction"` |

Pyrenees: `[tsgen] method = "series"`, `series_file =
"input/alpha_orbital_t2018.dat"`, `sigma = 0`; `[tsforcing] active = True`,
`f_ta = 10`, `f_to = f_so = 0`. tsgen samples the series at absolute time and
holds its end values; works today. Note that `tsforcing_init` runs tsgen even
when inactive, so the series file must exist.

## Work needed

### 1. tsgen normalisation (fesm-utils)

Done: on fesm-utils origin/dev 021f403 (the local clone checked was behind).

`norm_t0`/`norm_t1` are set in `yelmox_Greenland_snapesm*.nml` and the
CHANGELOG (Unreleased) says snapesm needs them, but they are not in tsgen
(nor in its history): a `series` index is the raw proxy (alpha_orbital_r2020:
1.04 at 0 ka, -13.8 at 21 ka), not 0 -> 1 between the snapshots.

- Add to `method = "series"`: `f = (s(t) - s(t0)) / (s(t1) - s(t0))` when
  both keys are set (no change otherwise).
- Needed by every former snap_1ind config (GRL, LIS, North, ANT paleo_lgp).

### 2. snapesm ocean

Done (forcing design phase 3): annual or monthly ocean, `so` optional.

`reduce_ocean` builds an ocean only when both `to` and `so` are bound, and
assumes a monthly source (`range`, `rep = 12`).

- Annual / time-less ocean (ANT reference `OCEAN_ISMIP6_J20`, Montoya
  `*-ocean`).
- Ocean without salinity (`so` = 0, as snapclim).
- Without these the ANT configs have no `now%to_ann` and `step_shelf` fails.

### 3. OBM air-temperature anomaly

`climate_air_anom` (`libs/yelmox_climate.f90`) works for snapclim only:
the raw `fname_at` series × `dTa_const`, at any time, independent of the
climate update. Options:

- (a) Backend-agnostic: the domain mean of `now%tsl_ann - ref%tsl_ann` from
  the last climate update. Works for snapesm, esm and rembo; snapesm's
  unused `dTa_const`/`dTo_const`/`dSo_const` are removed. Lags by `dt_clim`.
- (b) `snapesm_air_anom(sc, time)`: the `at` index × `dTa_const`, 0 without
  an `at` index. Closest to the current behaviour.

Proposed: (a). The rest of the OBM coupling (`obm2ism` writes
`dom%clim%now%to_ann`) is already backend-generic.

### 4. rembo ocean

Done (forcing design phase 3), with optional fields: the `pd` snapshot binds
only `to`/`so`, and the rembo backend adds `dT_ocn` itself.

rembo uses snapclim for its ocean (`const` + `dT_ocn`). Replace with
snapesm, `pd` snapshot only, `dTo = dT_ocn` passed to `snapesm_update`.
Caveat: snapesm always builds the atmosphere, so the `pd` snapshot needs
`zs`, `tas`, `pr` although rembo ignores them.

### 5. var_defs databases and Montoya from global files

Montoya2008 is in FesmData as `Paleoclimate-CLIMBER3a-<sim>` and
`-<sim>-ocean` (`present`, `lgm_1p7weak`, `lgm_1p7strong`), regional
products at 32 km only, global on GLOBAL-0.5/1/2/5DEG (ocean 1/2/5DEG).
Fields: `z_srf`, `t2m_ann`, `t2m_jja`, `pr_ann`; ocean `to(depth)`,
`mask_ocn`; units from the files. The remap already smooths (atmosphere
sigma 417 km, ocean 208 km; v1 files 250/100 km). On ANT-32KM the new
fields match the v1 ones closely (lgm_1p7strong t2m_ann mean -18.86 vs
-18.82 degC).

Read the global files with varslice online remap (already supported, incl.
lon-lat sources, 3D levels, map cache in `maps/`):

```fortran
&clim_c3a_lgm_strong_tas_ann
    filename = "fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-CLIMBER3a-lgm_1p7strong.nc"
    name     = "t2m_ann"
    units    = "degC" "K"
    scaling  = 1.0 273.15
    time     = 0 0 0 0 1
    remap    = "bilin"
    grid_src = "GLOBAL-1DEG"
/
&clim_c3a_lgm_strong_to
    filename = "fesmdata/Global/GLOBAL-1DEG/GLOBAL-1DEG_Paleoclimate-CLIMBER3a-lgm_1p7strong-ocean.nc"
    name     = "to"
    ...
    remap    = "bilin"
    grid_src = "GLOBAL-1DEG"
    fill     = "nearest"
/
```

Layout: one database per domain (`input/greenland_clim.nml` exists; new
`antarctica_clim.nml`, `north_clim.nml` for LIS + North, `pyrenees_clim.nml`,
`srg_clim.nml`), each with the domain's `pd` reference plus the generic PMIP3
(`{domain}/{grid_name}`) and global Montoya groups. Configs list only the
snapshots they use; the Montoya groups are available without being loaded.
Alternative: `var_defs` as a list of files, with one shared
`input/paleo_clim.nml` (no copies).

### 6. Dropped snapclim features

- the 80 km smoothing of clim0 before -1 ka (`f_corr`, snap_1ind/_new)
- `f_stdev` (precipitation variability), `f_p_ne`
- `f_hol` (already removed from snapesm)
- `hybrid`, `recon`, `miocene`, `*_abs`

### 7. Clean-up

- `libs/snapclim.f90`, its Makefile rules (`config/Makefile_yelmox.mk`),
  the snapclim cases in `libs/yelmox_climate.f90`.
- Stale tests: `tests/test_snapclim_ref.f90`, `tests/build_snapclim_ref.sh`
  (already broken: missing `dx`); `tests/test_snapesm_ref.f90` and
  `input/greenland_snp.nml` point at a removed worktree and the old schema.
- Dead code: `libs/obm/atm2ocean.f90` (`update_theta_from_snapclim`, never
  called).
- `yelmox/run_mg_forward.sh` overrides `snap.dTa_const`/`dTo_const`: move
  to `tsforcing.f_ta`/`f_to`. (Its par file has `fname_at = zeros.dat`, so
  the warm runs get 0 K today: existing bug.)
- Docs: `climate-snap.md` (snapclim section), `yelmox.md`,
  `yelmox-bipolar.md`, `climate-rembo.md`, `experiment-ant-paleo.md`,
  `multigrid.md`, `_quarto.yml`; CHANGELOG.
- `[comps] climate` default and the `snapclim` choice in the error messages
  (`kryos.f90`, `yelmox_climate.f90`).

### snapesm details to keep in mind

- Every key is required (`nml_read`): top level, per field (`kind`,
  `blend`, `index`), per snapshot (`monthly` and one entry per field).
- `zs`, `tas`, `pr` must be bound in every snapshot.
- `snapesm_end` is a stub; `dx` in `snapesm_update` is unused.
- Restart: one file per index; none with zero indices.

## Commits

1. fesm-utils dev: tsgen series normalisation (`norm_t0`, `norm_t1`).
2. yelmox: snapesm ocean (annual, no salinity); `climate_air_anom`; rembo on
   snapesm.
3. yelmox: var_defs databases with PMIP3 and global Montoya.
4. yelmox: config conversions, one commit per family: GRL (+ snapclim, pd,
   rembo), ANT (+ pd, paleo_spinup, paleo_lgp), LIS/North, Pyrenees
   (tsforcing), SRG, bipolar, `run_mg_forward.sh`.
5. yelmox: remove snapclim, tests, dead code; docs and CHANGELOG.

Then 10-yr Levante runs of each converted config.

## Open questions

1. OBM air anomaly: (a) from the climate output or (b) snapesm `at` index?
2. Databases: one per domain (copied PMIP3/Montoya groups) or a `var_defs`
   list?
3. Global resolution for Montoya: GLOBAL-1DEG (atmosphere and ocean) or
   finer for the atmosphere (GLOBAL-0.5DEG)? Bilinear or conservative?
4. dev only, or also `2.4-beta`?
