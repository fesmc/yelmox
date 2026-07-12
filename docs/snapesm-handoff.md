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

1. **fesm-utils fixes** (in progress / as needed):
   - tsgen `series_load_ascii`: add a plain-header skip (follow-up #2) so production `.dat` files
     load without the `tests/` fixture.
   - varslice: allocate the `ndim=4` missing-values branch in `varslice_update` (follow-up #3).
2. **Rebuild `include-serial`** from fesm-utils `dev` (its `varslice.mod` is current but `tsgen.mod`
   is **stale** — no `group` arg). ⚠ Do NOT rebuild while yelmox experiments run (project rule).
3. **Wire `climate_backend`** in `libs/yelmox_domain.f90` (mirror how `yelmox_esm` substitutes
   `esm_forcing_class`); add `snapesm.o` to the driver object lists in `config/Makefile_yelmox.mk`.
4. Validate in-model against snapclim, flip the default, retire `libs/snapclim.f90`.

Caller touch-points (from the design note): `yelmox_rembo` mutates `now%to_ann` and runs ocean-only;
`yelmox_bipolar/obm_coupling.f90` reaches into snapclim internals (`snp%at%time/%var`,
`par%dTa_const`) — needs index accessors + a `dTa_const` carrier (already in `par`).

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
