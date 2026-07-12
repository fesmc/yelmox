# snapesm — session handoff

Working handoff for the `snapclim` → `snapesm` refactor. Read with
[`snapesm-design.md`](snapesm-design.md) (design rationale + config model).

**Goal:** replace the monolithic `snapclim` (`libs/snapclim.f90`, ~2836 lines) with a thin
`snapesm` wrapper over the fesm-utils primitives `varslice` (field loading) and `tsgen`/`series`
(scalar forcing), validated to tight tolerance against snapclim on real configs.

**Status:** the numeric physics port is **complete and validated to tight tolerance** on
Greenland GRL-16KM `snap_1ind_new` + `fraction`. Remaining work is **in-model integration**
(wire a `climate_backend` switch into `yelmox_domain`, then retire snapclim).

---

## Where everything lives

**yelmox worktree** — `.claude/worktrees/snapclim2` (branch `snapclim2`, off `dev`; module `snapesm`).
- `libs/snapesm.f90` — the module (full physics; compiles + validates).
- `config/Makefile_yelmox.mk` — has the `snapesm.o` compile rule (not yet in any driver object list).
- `input/greenland_clim.nml` — var_defs database (per-file varslice metadata, pd/picontrol/lgm).
- `input/greenland_snp.nml` — assembly config (`&snp` + field/snapshot/index groups).
- `tests/test_snapesm_ref.f90`, `tests/test_snapclim_ref.f90` — the two validation harnesses
  (identical inline dumper, same synthetic z_srf/basins/times), `tests/diff_nc.py` — field diff,
  `tests/alpha_orbital_r2020.dat` — header-commented index fixture (see follow-up #2).

**fesm-utils** — symlink `fesm-utils -> /Users/alrobi001/models/fesm-utils`, checkout on **`dev`**
at `10ab52d`. The forcing API (`varslice subs`, `tsgen group`) is **merged to dev** (branch
`snp2-forcing-api` is now redundant). `dev` also carries the breaking condensed varslice namelist
(`units`/`scaling`/`time`) — commit `d35664d`.

**Reference outputs** — `logs/snapclim_ref.nc`, `logs/snapesm_ref.nc` (gitignored; regenerate below).

### Commits (branch `snapclim2`, newest first)
`d2523dc` validation test-bed · `e412632` numeric physics · `f5b87be` handoff · `9178566` doc→var_defs
· `ba4af91` rename→snapesm + var_defs model · `a71dd4e` reference harness · `326c25c` pipeline
structure · `9efdfe1` load via subs · `99260a1` par_load · `c48cbd9` skeleton.

---

## What is DONE

- **Design** (see design note): var_defs database + assembly config; named-member output state;
  generalized `snap(:)`/`idx(:)` with a distinguished `%ref`. **Revised during the port** to a
  per-field blend model: temp is linear (`ref + aa·(s1−s0)`), **precip is a ratio blend**
  (`ref·(aa·(s1/s0−1)+1)`), each field carries its own driving `index`, and snapshots carry a
  `monthly` flag (raw monthly vs annual+summer synthesis) and a `time` (for index normalization).
- **fesm-utils forcing API merged to `dev`** (varslice `subs`, tsgen `group`); condensed varslice
  namelist is the dev format.
- **snapesm physics (validated):**
  - `reduce_snapshot` (+ `reduce_ocean`, `model_depth`) — port of `read_climate_snapshot` +
    `read_ocean_snapshot`: monthly vs annual(+summer cosine synth) atm, sea-level lapse reduction,
    `prcor`, and the 42→23 ocean vinterp. Ocean is read with `rep=12` (whole-year sub-annual
    select) to avoid the sp/dp fractional-boundary trap.
  - `combine` (+ `blend_field`, `find_endpoints`, `norm_index`) — per-field linear/ratio blend;
    index normalized to snapshot times via `series_interp1`; Holocene `f_hol` scaling.
  - `transform` — tsl→tas inflation (phase `m*30-15`) + aggregates.
  - `derive` — `fraction` ocean, using the **pre-transform `ta_ann`** to reproduce snapclim's
    one-step lag exactly (decision: replicate faithfully now; fix as a separate commit later).
- **Validation** (`diff_nc.py` snapesm vs snapclim, GRL-16KM `snap_1ind_new`, times `0/-21000/-120000`):
  `tas,tsl,ta_ann,ta_sum,z_srf,mask` **bit-exact**; `pr,pr_ann,so_ann` ~1e-6; `to_ann` matches at
  all ocean cells except one deep level (800 m ≈ a coincident input node) at 0.73 K (accepted).

---

## What REMAINS — in-model integration

**Integration BUILDS end-to-end for BOTH backends.** All merged to `dev` (yelmox + fesm-utils);
`include-serial` rebuilt; Yelmo/FastIsostasy/yelmox libs rebuilt against current fesm-utils.
Verified: `make yelmox openmp=0` (snapclim, default) and `make yelmox CLIMATE=snapesm openmp=0`
both compile + link to `yelmox.x` (EXIT=0). **Remaining = run the model both ways and diff outputs.**

Done on branch `snapclim2` (now merged to `dev`):

- **fesm-utils fixes** (branch `snp-forcing-fixes`, off dev): tsgen/series plain-header skip;
  varslice `ndim=4` missing branch. Committed, not merged to dev.
- **Backend-agnostic climate** (`e412632`..`9febf81`):
  - `libs/climate_out.f90` — `climate_out_class {now, ref}`, each a `clim_state_class`
    {tas, pr, tsl_ann, ta_ann, pr_ann, to_ann, so_ann, depth}. Anomalies = now − ref.
  - `libs/yelmox_climate_snapclim.f90` / `_snapesm.f90` — same `module yelmox_climate` + API
    (`climate_init`/`climate_update`, fills `climate_out`); compile-verified vs each backend.
  - `libs/yelmox_domain.f90` — `ice_domain%snp` → `%cl` (adapter) + `%clim` (agnostic output);
    init/update via the adapter; all reads go through `dom%clim%now`/`%ref`.
  - `yelmox_rembo.f90`, `yelmox_bipolar/obm_coupling.f90` — updated to `dom%cl%snp` (backend
    internals; snapclim only) + `dom%clim%now%to_ann` (ocean writes).
  - `config/Makefile_yelmox.mk` — `CLIMATE ?= snapclim` selects adapter source + backend obj.

**Next session — run in-model & diff:**
1. Run a short yelmox simulation with the default (snapclim) `yelmox.x`, then rebuild with
   `make yelmox CLIMATE=snapesm` and run the SAME config. Author a `snp`/var_defs config for the
   in-model run (the offline `input/greenland_snp.nml` + `greenland_clim.nml` are the template;
   the in-model `&snap` group is `snap`+domain-suffix — the adapter passes `group="snap"//sfx`).
2. Diff the two runs' `snap.nc` / marine-shelf / smb outputs. Physics matches offline (atmosphere
   bit-exact; ocean matches at real cells), so any divergence is an integration bug.
3. Flip the default to snapesm, keep snapclim selectable a cycle, then retire `libs/snapclim.f90`
   + `yelmox_climate_snapclim.f90`.

Notes / caveats:
- The bipolar (`obm_coupling`) + `rembo` drivers reach into `dom%cl%snp` (snapclim internals), so
  they build only with `CLIMATE=snapclim`. `make yelmox_bipolar` / `yelmox_rembo` under
  `CLIMATE=snapesm` will fail until those are migrated off backend internals.
- Rebuilding `include-serial` to current fesm-utils `dev` (removed `ncio` `actual_range`) required a
  clean rebuild of Yelmo + FastIsostasy + the yelmox libs (their cached mods were stale). Already
  done this session; future fesm-utils bumps need the same.

---

## Build & run commands

The fesm-utils API is on `dev`, but `include-serial`'s `tsgen.mod` is stale, so snapesm still
builds against **scratch** varslice+tsgen (from `fesm-utils/src`) until `include-serial` is rebuilt.
Run from the yelmox **main tree** root so `ice_data/`, `input/`, `logs/`, `fesm-utils/` resolve.
Note **`wp = sp`** in this build; zsh does not word-split unquoted vars — use `${=VAR}`.

```bash
cd /Users/alrobi001/models/yelmox
S=$(mktemp -d)
INC=fesm-utils/include-serial          # current except tsgen.mod
FW=fesm-utils/src                       # dev source (has tsgen group + condensed varslice)
WT=.claude/worktrees/snapclim2
FF="-ffree-line-length-0 -O2 -g"        # NOTE: -0 not -none; -O0 -fcheck=all for debugging
NCF=$(ls -d /opt/homebrew/Cellar/netcdf-fortran/*/lib | head -1)
NCC=$(ls -d /opt/homebrew/Cellar/netcdf/*/lib | head -1)
NCLINK="-L$NCF -lnetcdff -L$NCC -lnetcdf -Wl,-rpath,$NCF -Wl,-rpath,$NCC"

# scratch varslice+tsgen from dev src, then snapesm, then link each harness
gfortran ${=FF} -I$INC -c $FW/varslice.f90 -o $S/varslice.o -J$S
gfortran ${=FF} -I$INC -c $FW/tsgen.f90    -o $S/tsgen.o    -J$S
gfortran ${=FF} -I$S -I$INC -c $WT/libs/snapesm.f90 -o $S/snapesm.o -J$S

gfortran ${=FF} -I$S -I$INC $WT/tests/test_snapesm_ref.f90 \
    $S/snapesm.o $S/varslice.o $S/tsgen.o -L$INC -lfesmutils ${=NCLINK} -o $S/test_snapesm_ref.x
$S/test_snapesm_ref.x                                   # -> logs/snapesm_ref.nc

gfortran ${=FF} -I$INC -c libs/snapclim.f90 -o $S/snapclim.o -J$S
gfortran ${=FF} -I$S -I$INC $WT/tests/test_snapclim_ref.f90 $S/snapclim.o -L$INC -lfesmutils ${=NCLINK} -o $S/test_snapclim_ref.x
$S/test_snapclim_ref.x                                  # -> logs/snapclim_ref.nc

python3 $WT/tests/diff_nc.py logs/snapesm_ref.nc logs/snapclim_ref.nc --rtol 1e-4 --atol 1e-3
```

---

## Follow-ups / known issues

1. **nml errors on any missing param** (`ERROR_NO_PARAM=.TRUE.`, not public, no setter). snapesm
   reads only the keys the config supplies (field: `kind/blend/index`; snapshot: `monthly/time/
   idx_coord`; `is_ref` derives from `ref_name`); `combine/enabled/*_flags` are struct-default only.
   **Decision: leave nml as-is for now.** A public toggle/query would restore that config surface.
2. **tsgen `series_load_ascii` rejects plain headers** — production `.dat` files have a bare
   `time value` header; the loader skips only `#`/`!` lines. snapclim's own `read_series` tolerated
   it. **Fix planned: add a header-skip in tsgen** (then drop the `tests/` fixture and point the
   config back at `input/alpha_orbital_r2020.dat`).
3. **varslice `ndim=4` missing branch** — `varslice_update`'s "indices not found" else-branch only
   allocates `vs%var` for ndim 1–3, so a 4-D field with no matching time slice segfaults on
   `vs%var = mv`. Worked around here (whole-year read always matches); **fix planned in varslice.**
4. **`to_ann` land-cell 273.15** — varslice guards missing (keeps `-9999`); snapclim adds K to the
   sentinel (`-9725.85`). Accepted (land cells); diff masks non-positive cells. If exact land match
   is ever needed, read ocean in °C and add 273.15 unconditionally in `reduce_ocean`.
5. **fraction-ocean lag** — snapclim reads `now%ta_ann` before recomputing it (block tagged
   `ajr: untested!!`). Replicated faithfully (decision (a)); the proper fix (use current-step
   `ta_ann`) is a deferred separate commit.
6. **800 m ocean interp edge** — 0.73 K at the one model level nearly coincident with an input
   depth node (799.55 m); `interp_linear` sp tie-break. Accepted as negligible.

## Cautions / discipline
- Worktree bash discipline: `git -C <worktree> …`; absolute paths for Read/Edit/Write.
- Don't rebuild `include-serial` (or merge fesm-utils→dev again) while yelmox experiments run.
- Tight-tolerance bar: watch the two cosine phases (`m*30-30` reduce vs `m*30-15` transform) and the
  south/north lapse sign; `south = (domain=="Antarctica")`.

## Open design items
- Field-registry set to name explicitly vs `extra(:)`; driver `dTa/dTo/dSo` vs internal `idx(:)`
  composition; restart provenance-record format (all deferred to integration/restart work).
