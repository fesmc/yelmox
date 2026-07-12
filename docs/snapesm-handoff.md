# snapesm — session handoff

Working handoff for continuing the `snapclim` → `snapesm` refactor. Read this together with
[`snapesm-design.md`](snapesm-design.md) (the design rationale + config model).

**Goal:** replace the monolithic `snapclim` (`libs/snapclim.f90`, ~2836 lines) with a thin
`snapesm` wrapper over the fesm-utils primitives `varslice` (field loading) and `tsgen`/`series`
(scalar forcing), validated to tight tolerance against snapclim on real configs.

---

## Where everything lives

**yelmox worktree** — `/Users/alrobi001/models/yelmox/.claude/worktrees/snapclim2`
(git worktree dir is named `snapclim2`; **branch `snapclim2`**, off `dev`; module is `snapesm`).
- `libs/snapesm.f90` — the new module.
- `config/Makefile_yelmox.mk` — has the `snapesm.o` compile rule (not yet in any driver's object list).
- `docs/snapesm-design.md`, `docs/snapesm-handoff.md` (this file).
- `tests/test_snapclim_ref.f90`, `tests/build_snapclim_ref.sh` — the validation-reference dumper.

**fesm-utils worktree** — `/Users/alrobi001/models/fesm-utils/.claude/worktrees/snp2-forcing-api`
(**branch `snp2-forcing-api`**, off `dev`). Two non-breaking API additions, committed but **NOT
merged to `dev`** and the shared `include-serial` build is **NOT rebuilt**:
- `src/varslice.f90` — `varslice_init_nml`/`varslice_par_load`/`parse_path` gained optional
  `subs(:,2)` for extra `{key}`→value path substitutions.
- `src/tsgen.f90` — `tsgen_init`'s `label` arg → `group`, used verbatim (default `"tsgen"`, no
  forced `tsgen_` prefix).
- `CHANGELOG.md` updated.

**Reference output** — `/Users/alrobi001/models/yelmox/logs/snapclim_ref.nc` (gitignored; regenerate
with the script below).

**Scratch build dir** (ephemeral, gone next session — must rebuild):
`<scratchpad>/fesmmod/` held `varslice.mod`, `tsgen.mod`, `snapesm.o`.

### Commits
Branch `snapclim2` (newest first): `9178566` doc→var_defs · `ba4af91` rename→snapesm + var_defs
model · `a71dd4e` reference harness · `326c25c` update-pipeline structure · `9efdfe1` load via
subs · `99260a1` par_load · `c48cbd9` skeleton.
Branch `snp2-forcing-api`: `10ab52d` forcing API (varslice subs + tsgen group).

---

## What is DONE

- **Exploration + agreed design** (see design note): var_defs database + assembly config;
  Option-A named-member output state (`allocated()` = in-use); generalized `snap(:)`/`idx(:)`
  with a distinguished `%ref`; five-knob unified model where `now = ref + Σ wₐ·snapshot` and each
  method is a **weight vector** (`snap_1ind` → `[0,-aa,+aa]`).
- **fesm-utils API additions** — `varslice subs` + `tsgen group` (branch `snp2-forcing-api`),
  compile-verified. Non-breaking (stock `varslice_init_nml` has no yelmox callers; `tsgen_init`
  has one no-arg caller in `yelmox_domain.f90`).
- **snapesm (compiles clean):**
  - Types: `field_spec_class`, `snapshot_spec_class`, `field_binding_class` (1–2 varslice
    `src(:)` per field), `snapshot_class` (`spec` + `bind(:)` + reduced `state`),
    `climate_state_class` (named members), `snapesm_param_class` (incl. `var_defs`), `snapesm_class`.
  - `snapesm_par_load` — top-level (`var_defs`, `combine`, `manifold`, `ref_name`, `lapse`,
    `dTa/dTo/dSo_const`), field registry (`&snp_field_<f>`), snapshot specs (`&snp_<snap>`),
    index names.
  - `snapesm_init(sc, filename, domain, grid_name, nx, ny, time, basins, group)` — note the
    added `time` arg (tsgen indices need it). Calls par_load → tsgen_init per index
    (`group=<base>_idx_<name>`) → load_snapshots.
  - `snapesm_load_snapshots` — per (snapshot, field) reads group ref(s) from `&snp_<snap>`
    (key = field name; 1 = monthly, 2 = `[ann,sum]`), loads each from the `var_defs` file via
    stock `varslice_init_nml`, primes (static once / time-varying `extrap`), then `reduce`; then
    `set_ref`.
  - `snapesm_update(sc, z_srf, time, domain, dTa,dTo,dSo, dx, basins)` — 5 stages:
    `advance_indices` (real), `refresh_loads` (real — re-slices + re-reduces time-varying
    sources), `weights`+`combine` (combine real; weights stub), `transform` (stub), `derive` (stub).
  - **Generic combine is real**: `accum2`/`accum3` do `now = ref + Σ w·snap` over allocated
    reduced-state members; `set_ref` copies the `is_ref` snapshot's state into `sc%ref`.
  - Stubs for `write_init`/`write_step`/`restart_write`/`restart_read`/`end`.
- **Validation reference** — `test_snapclim_ref.f90` runs snapclim on **Greenland GRL-16KM,
  `snap_1ind_new`** with synthetic-but-deterministic `z_srf`(dome)/`basins`(1.0) at times
  `[0, -21000, -120000]`, dumping `now%{...}` to `logs/snapclim_ref.nc`. Sanity numbers:
  `ta_ann`(mid) = 251.6 / 243.1 / 254.0 K (present / LGM / Eemian) — physically correct.

---

## What REMAINS (the numeric physics port, in order)

Everything below is stubbed in `snapesm.f90` with exact snapclim references in the stub comments.

### 1. `snapesm_reduce_snapshot(sc, s)` — the foundation
Extract each field's varslice `%var` into `snap(s)%state` (monthly), then reduce to sea level.
Also: annual→monthly synthesis (when `bind%src` has 2 entries = `[ann, sum]`), ocean vinterp,
`beta_p`.
- **BLOCKER first:** confirm the varslice `%var` layout empirically. Files:
  clim0 `tas(month,yc,xc)` monthly; clim1 `t2m_ann(yc,xc)` annual; ocn `to(month,depth=42,yc,xc)`.
  varslice stores `%var(nx,ny,nz,nt)`. Need to know where months/levels land (dim 3 vs dim 4) for
  a static vs time field. Read `varslice_update`/`varslice_init_data` extraction logic
  (`fesm-utils/src/varslice.f90`), or write a 10-line probe that loads one file and prints
  `shape(vs%var)`.
- **Reduction formulas** (snapclim.f90:1918-1948):
  - `tsl_ann = ta_ann + lapse(1)*z_srf`; `tsl_sum = ta_sum + lapse(2)*z_srf`
  - `pr = pr*(1 + f_stdev*pr_stdev_frac)`; `prcor_ann = pr_ann/exp(beta_p*(ta_ann-tsl_ann))`
  - monthly `tsl(m) = tas(m) + z_srf*(lapse(1) ± (lapse(2)-lapse(1))*cos(2π(m*30-30)/360))`
    (south `+(lapse2-lapse1)`, north `+(lapse1-lapse2)`; note phase **m*30-30** here)
  - `prcor(m) = pr(m)/exp(beta_p*(tas(m)-tsl(m)))`
- **Annual→monthly synthesis** (snapclim `clim_monthly=False` branch, ~1894-1910): cosine seasonal
  cycle from `ta_ann`+`ta_sum` (and `pr_ann`). Reproduce exactly.
- **`beta_p`** computed in `read_climate_snapshot` from `f_p`, `f_p_ne` (grep `beta_p`/`f_p` in
  snapclim; `f_p_ne` is hardcoded 1.0 for Greenland — snapclim.f90 ~1566/1579).
- **Ocean vinterp** (`read_ocean_snapshot`, snapclim.f90:1971-2168): read `to/so(month,depth,…)`,
  monthly-average, vertically interpolate 42-level input → the model's 23-level axis
  (nzo=23, 0–3000 m, snapclim.f90:300-311, `interp_linear`). Note the `depth0=abs(depth0)` fix.

### 2. `snapesm_weights(sc, w)` — port `snap_1ind_new`
- `calc_temp_1ind`: `temp_now = temp0 + aa*(temp2-temp1)` ⇒ `w = [0, -aa, +aa]` on
  `[ref=pd, s1=piControl, s2=lgm]` (ordered by `idx_coord`/registry).
- `snap_1ind_new` **index normalization** (snapclim.f90:437-459) — MUST read and port; it
  renormalizes the index `aa` to snapshot times. The index `at`/`ap` come from the tsgen series
  (`sc%idx`), value at `time` after `advance_indices`.
- Later methods (`snap_2ind`, `miocene`, `*_abs`, `hybrid`, `anom`, `const`, `recon`, `fraction`)
  are each another weight pattern / derive rule — port as needed per config.

### 3. `snapesm_transform(sc, z_srf)` — inflate + aggregates (snapclim.f90:779-812)
- `south = (domain=="Antarctica")`
- `tas(m) = tsl(m) - z_srf*(lapse(1) ± (lapse(2)-lapse(1))*cos(2π(m*30-15)/360))`
  (note phase **m*30-15** here, different from reduce's m*30-30)
- `pr(m) = prcor(m)*exp(beta_p*(tas(m)-tsl(m)))`
- `tsl_ann=mean(tsl)`, `ta_ann=mean(tas)`; `*_sum` over DJF `[12,1,2]` (south) / JJA `[6,7,8]`
  (north) /3; `pr_ann=mean(pr)/12*365`, `prcor_ann=mean(prcor)/12*365`.

### 4. `snapesm_derive(sc, dTa, dTo, dSo)`
- `fraction` ocean rule (ocean anomaly = `f_to`·mean atm anomaly) and folding the driver
  anomalies `dTa/dTo/dSo` into the fields (snapclim uses the arg when present, else its own index).

### 5. Config files (author these)
- `input/greenland_clim.nml` — the var_defs database: varslice groups for pd/piControl/lgm ×
  {tas(+sum), pr, zs, to, so}. See the design note for the exact key layout and the real file paths
  (all present under `ice_data/Greenland/GRL-16KM/`). Reference snapclim's `&snap_clim0..2` in
  `yelmox/yelmox_Greenland.nml` for which file/var each maps to.
- The `&snp…` assembly config — put it in a test nml (or extend a copy of the Greenland nml).
  `snapshots="pd" "picontrol" "lgm"`, `ref_name="pd"`, `fields="tas" "pr" "to" "so" "zs"`,
  `indices="at" "ap"`, `lapse` and `f_*` copied from `&snap` in `yelmox_Greenland.nml`.

### 6. snapesm harness + diff
- Copy `tests/test_snapclim_ref.f90` → `tests/test_snapesm_ref.f90`, swap `snapclim`→`snapesm`
  (add the `time` init arg), **same synthetic z_srf/basins and times**, write to
  `logs/snapesm_ref.nc`.
- Diff `logs/snapesm_ref.nc` vs `logs/snapclim_ref.nc` field-by-field (a small Python/`ncdump`
  or Fortran max-abs-diff) on `tas/pr/tsl_ann/ta_ann/pr_ann/to_ann/so_ann`. Iterate to tight
  tolerance. Start with `const`/`anom` sanity, then `snap_1ind_new`.

### 7. Integration (after validation)
- Merge `snp2-forcing-api` → fesm-utils `dev`, rebuild `include-serial` (see caution below).
- Wire a `climate_backend = "snapclim" | "snapesm"` switch in `libs/yelmox_domain.f90`
  (mirror how `yelmox_esm` substitutes `esm_forcing_class`); add snapesm.o to the driver object
  lists in `config/Makefile_yelmox.mk`. Validate in-model, then retire snapclim.

---

## Build & run commands

**All builds are manual against a scratch dir** because the fesm-utils API changes are not merged.
Run from the yelmox **main tree** root (`/Users/alrobi001/models/yelmox`) so relative paths
(`ice_data/`, `input/`, `logs/`, `fesm-utils/`) resolve — the worktree does not have the
gitignored `ice_data`/`input`/`fesm-utils` symlinks.

```bash
cd /Users/alrobi001/models/yelmox
S=$(mktemp -d)                                   # scratch mod/obj dir
INC=fesm-utils/include-serial                    # built shared deps (precision/ncio/nml/series/constants/mapping)
FW=/Users/alrobi001/models/fesm-utils/.claude/worktrees/snp2-forcing-api/src
WT=.claude/worktrees/snapclim2
FF="-ffree-line-length-none -O2"
NCF=/opt/homebrew/Cellar/netcdf-fortran/4.6.2/lib
NCC=/opt/homebrew/Cellar/netcdf/4.10.0/lib
NCLINK="-L$NCF -lnetcdff -L$NCC -lnetcdf -Wl,-rpath,$NCF -Wl,-rpath,$NCC"

# (a) rebuild the modified fesm-utils modules into scratch (override the shared lib versions)
gfortran $FF -I$INC -c $FW/varslice.f90 -o $S/varslice.o -J$S
gfortran $FF -I$INC -c $FW/tsgen.f90    -o $S/tsgen.o    -J$S

# (b) compile snapesm against the NEW mods (scratch first, then shared)
gfortran $FF -I$S -I$INC -c $WT/libs/snapesm.f90 -o $S/snapesm.o -J$S

# (c) link a snapesm harness: put scratch varslice.o/tsgen.o BEFORE -lfesmutils so they win
gfortran $FF -I$S -I$INC $WT/tests/test_snapesm_ref.f90 \
    $S/snapesm.o $S/varslice.o $S/tsgen.o -L$INC -lfesmutils $NCLINK -o $S/test_snapesm_ref.x
$S/test_snapesm_ref.x

# --- the snapclim REFERENCE (already produces logs/snapclim_ref.nc) ---
gfortran $FF -I$INC -c libs/snapclim.f90 -o $S/snapclim.o -J$S
gfortran $FF -I$S -I$INC $WT/tests/test_snapclim_ref.f90 $S/snapclim.o -L$INC -lfesmutils $NCLINK -o $S/test_snapclim_ref.x
$S/test_snapclim_ref.x
```

(The committed `tests/build_snapclim_ref.sh` encodes the reference build but assumes cwd has both
`tests/` and the data; run the explicit commands above from the main tree if in doubt.)

---

## Cautions / discipline

- **Worktree bash discipline** (project rule): use `git -C <worktree> …`; use absolute/worktree
  paths for Read/Edit/Write; a bare `git status` in the main tree does NOT reflect the worktree.
- **fesm-utils is a shared production lib.** `snp2-forcing-api` is committed but **unmerged**;
  the shared `include-serial` is unchanged, so **snapclim still builds normally** and any running
  experiments are unaffected. snapesm only builds against the scratch mods (or after a merge+rebuild).
  Do NOT merge + rebuild `include-serial` while batch experiments run on it (project rule).
- Keep the fesm-utils branch name `snp2-forcing-api` (don't rename for the module rename).
- **Tight-tolerance bar** (user requirement): port the `calc_*` formulas exactly; validate every
  method against `snapclim_ref.nc` before trusting it. Watch the two different cosine phases
  (`m*30-30` in reduce vs `m*30-15` in transform) and the south/north lapse sign.
- Reproducibility of the reference does not need the fesm-utils changes (snapclim uses the old API).

## Open design items (from the design note)
- Final field-registry set to name explicitly vs `extra(:)`.
- Whether driver `dTa/dTo/dSo` compose with snapesm's own `idx(:)` or are mutually exclusive per field.
- Provenance-record format for restart (and whether to snapshot derived fields for verification).
