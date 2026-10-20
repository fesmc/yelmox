# Par files: bases plus -p, forcing library

Status: plan (2026-10-20), not started. The leading `_` keeps quarto from
rendering it.

States checked: yelmox dev 6ec65b4.

## Motivation

Survey findings:

- `yelmox/` holds 26 par files. Ten are bases: nine domain × backend files
  (Greenland, Antarctica, North, LIS, Pyrenees, SRG, esm_Greenland,
  esm_Antarctica, rembo_Greenland) plus GreenlandPaleo. The other 16 are 15
  variants of them and the rembo aux file `rembo_Greenland.nml`.
- Variants are of two kinds.
  - Flat variants differ only in key values: times, grid, `esm.par_file`,
    `surface_method`, data paths. These are esm_*_1pctCO2,
    esm_Greenland_tipmip, esm_Greenland_k24, esm_Antarctica_ismip7,
    esm_Antarctica_nudge, pd_Antarctica, Antarctica_paleo_spinup and
    GreenlandPaleo_spinup.
  - Structural variants add or remove named groups: `snap_<name>`,
    `snap_axis_<name>`, `idx_<name>`, `ytrc` plus `elsa_*`, `negis`, `bias`,
    `taper`. These are GreenlandPaleo_lgp, _lgm, _deglac,
    Antarctica_paleo_lgp, Greenland_ramp and pd_Greenland.
- Cost. In the 120 days before 2026-10-20, 148 of 541 commits touched par
  files. Each variant rode along on 45 to 50 of them, for renames and new
  keys. Several scripts (run_ismip7, run_tipmip, run_1pctco2) already re-set
  with `-p` what their variant file changed, so file and script duplicate
  each other.
- Named groups are duplicated across files. `snap_picontrol` and `snap_lgm`
  are identical in Greenland, GreenlandPaleo_lgp, North and LIS; Antarctica's
  differ. `idx_a`, `idx_at` and `snap_axis_a` appear in 4 files each.
  `snap_pd` is domain-specific.
- runme `-p` can set any key present in the par file, including vectors with
  `==` and ensembles. It cannot create a group: a key absent from the par
  file is an error unless a par_defaults file declares it, and
  `.runme/info.json` does not set par_defaults.
- The Fortran reads `snap_<name>`, `snap_axis_<name>` and `idx_<name>` from
  the par file only. Varslice groups are already searched through a list of
  files: `var_defs_find` in `libs/yelmox_nml.f90` looks in the var_defs files
  first, then the par file.

Considered and dropped:

- Case files, i.e. partial namelist overlays merged by runme. Once the
  library exists, their only unique capability, adding groups, is gone. They
  add a third layer, go stale like any partial namelist, and are not
  validated until merged. A `-p` key that was renamed fails loudly at runme
  time, which is better.
- A fixed par-file name, `yelmox.nml` in every run dir. Deferred: distinct
  names make it easier to trace a run's origin, and it is more code change
  than it is worth now.

## Principle

`yelmox/` holds one par file per domain × climate backend. Everything an
experiment changes is a `-p` list in its script. Shared named groups
(snapshots, axes, indices, elsa) live once in a library file.

Bases kept (10): Greenland, GreenlandPaleo, Antarctica, North, LIS,
Pyrenees, SRG, esm_Greenland, esm_Antarctica, rembo_Greenland. The rembo aux
file `rembo_Greenland.nml` stays too.

- GreenlandPaleo stays because a different grid hub and border treatment is
  a domain, not an experiment.
- esm_Antarctica keeps its current data setup: BedMachine, Zwally basins,
  UKESM. ISMIP7 is a variant expressed in `scripts/ismip7`.

Deleted (15):

- GreenlandPaleo_spinup, _lgp, _lgm, _deglac
- Greenland_ramp, pd_Greenland
- pd_Antarctica, Antarctica_paleo_spinup, Antarctica_paleo_lgp
- esm_Greenland_1pctCO2, _tipmip, _k24
- esm_Antarctica_1pctCO2, _ismip7, _nudge

## Fortran

- New library file `input/forcing_db.nml`:
  - `snap_<name>`: picontrol, lgm, the Antarctic PMIP3 pair with suffix
    `_ant`, b18_base, b18, b20.
  - `snap_axis_<name>`: a, at, b18, b20.
  - `idx_<name>`: a, at, one, ramp, orb, hyst, dT.
  - `elsa_<name>`: grl_paleo, ant_paleo.

  `snap_pd` stays in each base as the domain's reference. `[snap]`,
  `[snap_field_*]` and `[indices]` stay in the base, since they are the
  per-run selection.
- Lookup rule, one for all groups: the par file first, then the var_defs
  files in order. This flips the current rule in `var_defs_find`. A par file
  that redefines a library group wins.
- `libs/snapshots.f90`: `snap_<name>` and `snap_axis_<name>` are resolved
  through `var_defs_find` with the `[snap]` var_defs list, which gains
  `input/forcing_db.nml`.
- `libs/kryos_forcing.f90`: new `[indices]` var_defs key, default empty,
  declared in `input/yelmox_defaults.nml`. `idx_<name>` is resolved the same
  way.
- elsa: `ytrc.elsa_nml = "input/forcing_db.nml"`. No code change.
- `scripts/check_par_nml.py` follows the library for named groups. It
  validates every group of `forcing_db.nml` against its template by prefix:
  `snap_axis_*` -> `snap_axis`, `snap_*` -> `snap_snapshot`, `idx_*` -> `idx`.

## runme and .runme/info.json

- par_defaults: `input/yelmo_defaults.nml` and `input/yelmox_defaults.nml`,
  so `-p bias.active=True` works on a base without `[bias]`.
- New runme feature `par_library`, a list of namelist files in `info.json`.
  When a `-p` key names a group absent from the run's par files but present
  in a library file, runme first copies the whole group into the staged par
  file, then applies the override. This keeps `-p idx_a.f_min=0` working
  after `idx_a` moves to the library, and the staged file stays
  self-contained for that group. The Fortran rule above then needs no
  key-level merging.
- Precedence for a `-p` key:
  1. group in the par file (normal);
  2. else the whole group from `par_library`;
  3. else the single key from par_defaults;
  4. else error.
- Nothing else changes: `-n` is the base, staged under its own name and
  passed as the argument.

## Scripts and docs

- Each deleted variant becomes a commented `-p` array in its script, in the
  style of `scripts/ismip7/opt_grl_ismip.sh`.
- Arrays used by more than one script go in a sourced `params.sh` in that
  experiment folder:
  - `scripts/grl-paleo`: the b18/b20 reconstruction forcing shared by the lgm
    spin-up and the deglaciation.
  - `scripts/ant-paleo`: the 32 km physics.
  - `scripts/ismip7`: the ISMIP7 Antarctic data sources, also sourced by
    `scripts/1pctCO2/opt_ant.sh`.
- Rule of thumb: experiment-defining settings go in the arrays, with a
  comment each. Run control goes in the invocation: run_step, times,
  restart, grid, output steps, queue.
- pd_Greenland, pd_Antarctica and Greenland_ramp become one-line commands in
  the experiment docs.
- esm_Antarctica_nudge has no script today and gets
  `scripts/nudge/run_nudge_antarctica.sh`.
- Case tables in `scripts/docs-figures/run.sh` and
  `scripts/climate-compare/run.sh`: the base in the par column, the variant's
  arrays in the overrides column.
- `scripts/run_chunks.sh` and `yelmox/run_mg_*.sh` are unchanged.
- Docs: experiment pages (Config lines), `docs/yelmox.md` Configuration,
  `docs/climate-snapshots.md` (lookup rule, library), CHANGELOG.

## Verification

- `check_par_nml` on every base and on the library.
- The 18 docs-figures cases on Levante. The staged par files must equal the
  deleted variants up to comments, so the climate output must match the
  current figures' data exactly. Compare the NetCDF fields, not the figures.
- Staged-only runs of the ISMIP7 ANT, 1pctCO2 and TIPMIP workflows, to
  inspect the staged files. Plus one short real spin-up step for ISMIP7 ANT.

## Order

1. runme: `par_library`, version bump, install on Levante.
2. yelmox Fortran: lookup rule, snapshot and index lookup, new key,
   `check_par_nml`.
3. `forcing_db.nml`, trimmed bases, `info.json`.
4. Scripts, deletions, case tables. Test on Levante.
5. Docs and CHANGELOG.
