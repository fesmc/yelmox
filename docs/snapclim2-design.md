# snapclim2 — design note for the snapclim replacement

Design specification for replacing the monolithic `snapclim` package (`libs/snapclim.f90`,
~2836 lines) with a thin orchestration layer built on the fesm-utils primitives `varslice`
(NetCDF field loading) and `tsgen`/`series` (scalar transient forcing). In progress.

> **Name.** `snapclim2` (module `snapclim2`, type `snapclim2_class`, file `libs/snapclim2.f90`)
> is a provisional placeholder so the new package can coexist with `snapclim` during migration.
> Rename before finalizing if preferred.

## Goal

Replace snapclim with a package that:

1. Delegates **all field loading** to `varslice` and **all time-series** to `tsgen`/`series`,
   deleting snapclim's hand-rolled NetCDF readers (`read_climate_snapshot`,
   `read_ocean_snapshot`, the two `*_reconstruction` variants) and series reader
   (`read_series`) — roughly the bulk of the 2836 lines.
2. Exposes **one clean, flexible API** in which today's ~11 forcing methods
   (`const/anom/snap_1ind[_new]/snap_1ind_miocene/snap_2ind/*_abs/hybrid/recon/fraction`)
   are *configurations* of a single blend model, not separate code paths — so new scenarios
   are namelist changes, not new kernels.
3. Handles **variable numbers of climate snapshots and driving indices** (was hardwired
   `clim0..clim3` + `at..bs`), with a **distinguished reference state**.
4. Handles fields that are **static or time-varying** with the same knob (including surface
   elevation `z_srf`, which some GCMs now vary in time).
5. Optionally loads **exotic forcing fields** (`smb`, `bmb_shlf`, …) as boundary forcing —
   off by default, no change to the standard climate path.
6. Supports **read/write and restart**, including a **provenance record** that documents the
   forcing state associated with a restart even when nothing needs to be reloaded.

## Non-goals

- Not changing what smbpal / marine_shelf compute. snapclim2 produces climate (and optionally
  exotic) *boundary forcing*; SMB/BMB stay in their existing modules unless a run explicitly
  loads them as forcing.
- Not a big-bang swap. Built alongside snapclim and switched in behind a flag (as
  `yelmox_esm` already substitutes `esm_forcing_class`), validated, then snapclim retired.

---

## Layered architecture

```
  Loading   varslice (fesm-utils)   one varslice_class per (snapshot × field)
  Indices   tsgen / series          one tsgen_class per driving index  (+ restart, kill)
  Wrapper   snapclim2 (NEW)         blend + transforms + derive + output state + restart/provenance
```

Only the wrapper is new code. Pattern precedents already in the tree:

- **`tsforcing_class`** (`libs/yelmox_domain.f90:154`) — the shape to copy: a small wrapper
  that owns the fesm-utils generator, reads its own namelist, exposes derived outputs, and
  delegates restart/kill to the library.
- **`varslice_init_nml_esm` / `parse_path_esm`** (`libs/esm.f90:855`) — filename templating
  (`{domain}`, `{grid_name}`, and extra placeholders) layered on stock varslice.
- What to **avoid** from esm: ~40 hand-named `varslice_class` members and a 3-way split
  update surface. snapclim2 uses arrays/registries and a single `update`.

---

## The unified model (five orthogonal knobs)

Every output field is:

> **`field = combine( snapshots, weights )`  → `transform`  → `derive`**

driven by five independent knobs. Every existing snapclim method is a setting of these:

| knob | values | absorbs |
|------|--------|---------|
| **snapshot kind** | `gridded-static` · `gridded-time` · `uniform` | recon (→ gridded-time), hybrid/anom (→ uniform) |
| **index manifold** | 0-D · 1-D · M-D, N snapshots | 1ind, 2ind, miocene (1-D, N=4) |
| **combine mode** | `anomaly` (vs ref) · `absolute` | `*_abs` variants |
| **derive rules** | e.g. `ocn = f·atm_anom` | `fraction` (ocean) |
| **transforms** | lapse · precip-scaling · seasonal-synthesis · unit (unit via varslice) | the `calc_*` kernels |

Method → configuration mapping:

| snapclim mode | snapclim2 configuration |
|---|---|
| `const` | one snapshot = ref, weight 1 |
| `anom` | uniform snapshot, offset from an index; anomaly combine |
| `snap_1ind`, `_new` | N gridded snapshots, 1-D manifold, anomaly combine |
| `snap_1ind_miocene` | same, N=4 |
| `snap_2ind` | gridded snapshots on a 2-D manifold |
| `snap_1ind_abs`, `snap_2ind_abs` | as above, **absolute** combine |
| `hybrid` | uniform snapshot carrying a 12-channel monthly series (tsgen `nc=12`) |
| `recon` | a **gridded-time** snapshot (varslice `with_time` + `interp`); anomaly combine — deletes the reconstruction readers |
| `fraction` (ocn) | a post-blend **derive** rule |

The two ~130-line `select case` dispatch blocks and the ~20 near-duplicate `elemental calc_*`
kernels collapse into: a generic weighted combine, a small set of shared transform helpers,
and an optional derive step.

---

## Types

Sketch (final field names/kinds settle during implementation; single precision preserved).

```fortran
! ---- configuration (immutable after load) ----------------------------------

type field_spec_class                 ! one physical field in the registry
    character(len=32)  :: name        ! "tas","pr","to","so","zs","smb","bmb_shlf",...
    character(len=16)  :: kind        ! "atm_monthly" | "ocn_3d" | "scalar_2d"
    character(len=16)  :: combine     ! "anomaly" | "absolute"
    logical            :: enabled     ! off by default for exotic fields
    ! transform flags: apply_lapse, precip_scaling, seasonal_synth, ...
end type

type snapshot_spec_class              ! one climate state: PD, LGM, piControl, ...
    character(len=64)  :: name
    logical            :: is_ref      ! the distinguished reference state
    real(wp)           :: idx_coord(2)! position on the 1-D/2-D index manifold
    ! per-field varslice configs are built from field_spec + this snapshot's paths
end type

! ---- loaded input ----------------------------------------------------------

type snapshot_class
    type(snapshot_spec_class) :: spec
    type(varslice_class), allocatable :: fld(:)   ! one per enabled field (registry order)
end type

! ---- output state (Option A: named members, unallocated when unused) -------

type climate_state_class
    ! atmosphere (monthly nx,ny,12)
    real(wp), allocatable :: tas(:,:,:), pr(:,:,:), sf(:,:,:), pr_stdev_frac(:,:,:)
    ! atmosphere annual/summer (nx,ny)
    real(wp), allocatable :: ta_ann(:,:), ta_sum(:,:), pr_ann(:,:)
    real(wp), allocatable :: tsl_ann(:,:), tsl_sum(:,:), prcor_ann(:,:)
    ! ocean
    real(wp), allocatable :: to_ann(:,:,:), so_ann(:,:,:)
    real(wp), allocatable :: depth(:)
    ! geometry / masks
    real(wp), allocatable :: mask(:,:), z_srf(:,:)
    ! exotic / expansion forcing (named for the common ones)
    real(wp), allocatable :: smb(:,:), bmb_shlf(:,:,:)
    ! escape hatch for genuinely ad-hoc future fields (temp-copy access only)
    type(named_field_class), allocatable :: extra(:)   ! %v("name") returns a copy
end type

! ---- main object -----------------------------------------------------------

type snapclim2_class
    type(snapclim2_param_class) :: par
    type(field_spec_class),    allocatable :: registry(:)  ! field definitions
    type(snapshot_class),      allocatable :: snap(:)      ! was clim0..clim3
    type(tsgen_class),         allocatable :: idx(:)       ! was at..bs (each may be nc-channel)
    type(climate_state_class)  :: now                      ! current output
    type(climate_state_class)  :: ref                      ! distinguished reference (was clim0)
end type
```

Notes:

- **Option A confirmed.** Named members are indexed directly (`now%tas(1,1)`); Fortran forbids
  subscripting a function result, so a generic `%v("tas")(1,1)` is impossible — the `extra(:)`
  map is only for rarely-touched fields where a temporary copy is acceptable. Unused fields
  are left **unallocated**; `allocated(now%smb)` is the "is this run using it" guard.
- **`%ref` is a real state**, not a view: callers hold references to it, and absolute modes need
  a well-defined baseline (replaces snapclim's `clim0 = clim1` trick).
- **Snapshot & index counts are variable.** `snap(:)` and `idx(:)` are allocatable; the fixed
  `clim0..3` / `at..bs` layout is gone.

---

## Namelist design

Base group name is passed via `group=` (default `snapclim2`, or `snapclim2_1`/`_2` for
multi-domain), mirroring snapclim's `group=`.

```
&snapclim2                     ! top-level
    combine     = "anomaly"    ! default combine mode
    manifold    = 1            ! index-manifold dimension (0/1/2)
    ref_name    = "pd"         ! which snapshot is the reference
    lapse       = ...
    ! precip/ocean scaling factors, f_stdev, etc.
    snapshots   = "pd" "lgm"          ! ordered list -> &snapclim2_snap_<name>
    fields      = "tas" "pr" "to" "so" "zs"   ! registry -> &snapclim2_field_<name>
    indices     = "at" "ap"           ! -> &snapclim2_idx_<name> (tsgen groups)
/

&snapclim2_field_tas
    kind        = "atm_monthly"
    combine     = "anomaly"    ! override top-level
    apply_lapse = True
/

&snapclim2_field_smb           ! exotic forcing, off unless enabled
    kind        = "scalar_2d"
    enabled     = False
/

&snapclim2_snap_pd             ! reference snapshot
    is_ref      = True
    idx_coord   = 0.0
    ! per-field file bindings (see below)
/

&snapclim2_snap_lgm
    idx_coord   = 1.0
/

&snapclim2_idx_at              ! a driving index == a tsgen group
    method      = "series"
    series_file = "..."
    series_var  = "..."
/
```

**Per-snapshot × per-field file binding.** Each snapshot needs a varslice config per enabled
field. To avoid snapclim's `clim_names(4)` positional array and esm's ~40 named members, bind
them as a small sub-group per (snapshot, field), e.g. `&snapclim2_snap_pd_tas`, reusing the
stock varslice namelist keys via a `varslice_init_nml_snapclim2` wrapper (adds any extra
`{...}` placeholders on top of `{domain}`/`{grid_name}`).

**varslice namelist convention (condensed form).** snapclim2 consumes whatever `varslice_init_nml`
accepts; the current condensed convention is:

```
&snapclim2_snap_pd_tas
    filename = "ice_data/{domain}/{grid_name}/.../..._monthly.nc"
    name     = "t2m"
    units    = "K" "K"                        ! units_in units_out
    scaling  = 1.0 0.0                         ! unit_scale unit_offset
    time     = 1.0, 1979.0, 2022.0, 1.0, 12.0 ! [active?, t0, t1, dt, n_sub]
/
```

The **first `time` value flags whether the time axis is active** (replaces the separate
`with_time` key); the remaining four are the old `time_par` `[t0, t1, dt, n_sub]`. A snapshot
with time-varying elevation just sets its `_zs` sub-group's `time` active flag to 1 — no code
path change. `units`/`scaling` collapse the former `units_in`/`units_out` and
`unit_scale`/`unit_offset` pairs onto one line each.

---

## Update algorithm

`update(sc, z_srf, time, domain, dTa, dTo, dSo, dx, basins)` — same signature as
`snapclim_update`, so the `yelmox_domain` call sites are unchanged.

```
1. Advance indices:  for each idx(k): tsgen_update(idx(k), time, var)   -> f_now (weights source)
                      (dTa/dTo/dSo, when passed, override/augment index-derived offsets)
2. Refresh loads:     for each snapshot s, enabled field f:
                          varslice_update(snap(s)%fld(f), time, method)   ! static -> no-op after 1st
                          ! gridded-time snapshots (recon, time-varying zs) slice/interp here
3. Combine:           for each enabled field f:
                          weights = manifold_weights(idx, snap%idx_coord)
                          if combine=="anomaly": now%f = ref%f + Σ w_s·(snap(s)%f - ref%f)
                          else (absolute):        now%f = Σ w_s· snap(s)%f
                          ! uniform snapshots contribute a spatially-constant term
4. Transform:         lapse-rate reduce/inflate (tsl<->tas using source & target z_srf),
                      precip-temperature scaling, seasonal-cycle synthesis for annual-only
                      snapshots, aggregates (ta_ann/ta_sum/pr_ann, ocean vertical interp to depth)
5. Derive:            optional rules, e.g. ocean anomaly = f_to · mean(atm anomaly)  [was `fraction`]
```

Steps 3–5 are shared helpers keyed by `field_spec%kind`, replacing the per-method kernels.
Static snapshots make step 2 a cached no-op after the first call (varslice already
short-circuits unchanged updates).

---

## Elevation (`z_srf`) handling

Two distinct elevations, previously conflated:

- **Source elevation** — where a snapshot's climate was defined; a normal snapshot field
  (`kind="scalar_2d"`, e.g. registry name `zs`). Time-inactive → constant (old behavior);
  time-active → time-varying (transient-orography GCM). The lapse transform reads the
  **time-current** source elevation each `update`, so `tsl` is recomputed when it varies.
- **Target elevation** — the model's current surface, the `z_srf` argument to `update`; used to
  inflate `tsl → tas`.

Only change vs today: the lapse-reduction helper reads `snap%zs` at update time instead of
assuming it is frozen at load. Localized to one helper.

---

## Output contract & caller compatibility

Minimum output set consumed by `yelmox_domain` (must be produced, same names):

- Atmosphere → smbpal: `now%tas`, `now%pr`; → smb_simple: `now%tsl_ann`.
- Ocean → marine_shelf: `now%to_ann`, `now%so_ann`, `now%depth`; anomaly baseline
  `ref%to_ann` (was `clim0%to_ann`).
- Diagnostics (`snap.nc`): `now%ta_ann`, `now%pr_ann`; `ref%ta_ann`.

Compatibility shims for the two special callers:

- **`yelmox_rembo`** mutates `now%to_ann` directly and calls update for the ocean only — keep
  `now%to_ann` writable and allow atmosphere-disabled configs.
- **`yelmox_bipolar/obm_coupling.f90`** reads `snp%at%time`, `snp%at%var`, `snp%par%dTa_const`
  and writes `now%to_ann`. Provide an index accessor exposing an index's `(time, var)` series
  (the `tsgen`'s underlying `series`), and carry a `dTa_const`-equivalent in `par`. Flag this
  as an explicit migration touch-point — obm is the one place reaching into snapclim internals.

---

## Read/write & restart

- **Diagnostics.** `snapclim2_write_init` / `snapclim2_write_step` write the `now` fields to a
  NetCDF file, generalized to the registry (loop over enabled fields) instead of the hardcoded
  list in `snapclim_write_step`.
- **Restart — prognostic state.** The only genuine prognostic state is the driving indices →
  delegate to `tsgen_restart_write` / `tsgen_restart_read` per `idx(:)` (handles Δf, kill,
  controller history). Snapshots are static input → not serialized; on read, indices are
  restored and all fields recomputed from `time`.
- **Restart — provenance record (the explicit ask).** Alongside the tsgen restart, write a
  self-contained record documenting the forcing state at the restart:
  - snapshot inventory: names, source files, `idx_coord`, `is_ref`;
  - enabled field registry + combine modes + transform flags;
  - index values and derived weights at restart time;
  - optionally a snapshot of the derived `now` fields for bit-tolerance verification.
  On read this is **documentation/validation only** — not reloaded to drive the run (fields are
  recomputed from the restored indices), satisfying "even just to document the state."
- Add `snapclim2_restart_write/read` to the `domain_restart_write/read` bundle in
  `yelmox_domain` (snapclim contributes nothing there today).

---

## Flexibility / expansion examples

- **New scenario, existing fields** — add a `&snapclim2_snap_<name>` group with its `idx_coord`
  and file bindings; extend `snapshots`. No code.
- **Higher-dimensional forcing** — set `manifold = 2` and give snapshots 2-D `idx_coord`; add a
  second index. No code.
- **Load SMB as forcing** — set `enabled = True` on `&snapclim2_field_smb`, bind files per
  snapshot; `now%smb` populates. Downstream can choose it over smbpal. No type change.
- **Time-varying orography** — set the time active flag in a snapshot's `_zs` sub-group.
- **New exotic field** — add to the registry; named member if hot, else `extra(:)`.

---

## Migration & validation plan

1. **Build alongside.** Add `libs/snapclim2.f90`; do not touch `snapclim.f90`. Wire into
   `yelmox_domain` behind a switch (e.g. `climate_backend = "snapclim" | "snapclim2"`),
   mirroring the `yelmox_esm` substitution pattern. Work in a git worktree (per repo policy).
2. **Port physics to tight tolerance.** Re-express the `calc_*` kernels as shared transform
   helpers; keep formulas numerically equivalent. Build a validation harness that runs a set of
   representative live configs (Antarctica/Greenland, at least one of each method in use) under
   both backends and diffs `now%{tas,pr,tsl_ann,to_ann,so_ann}` field-by-field against a tight
   tolerance. Treat any exceedance as a port bug.
3. **Cover each method in use.** Confirm which of `const/anom/snap_1ind[_new]/miocene/snap_2ind/
   *_abs/hybrid/recon/fraction` appear in active `.nml` files; validate each. Drop nothing until
   its replacement configuration passes.
4. **Restart parity.** Verify write→read→continue reproduces a no-restart run, and that the
   provenance record round-trips.
5. **Switch default & retire.** Flip the default to `snapclim2`, keep snapclim selectable for
   one cycle, then remove `libs/snapclim.f90` and its references.

## What gets deleted (once validated)

- `read_climate_snapshot`, `read_ocean_snapshot`, `read_climate_snapshot_reconstruction`,
  `read_ocean_snapshot_reconstruction` → `varslice`.
- `read_series`, `series_interp`, the `series_type`/`series_2D_type` handling → `tsgen`/`series`.
- The ~20 `elemental calc_*` kernels → a handful of shared transform helpers.
- The two `select case(atm_type)` / `select case(ocn_type)` dispatch blocks → generic combine.
- Fixed `clim0..clim3` + `at..bs` layout → `snap(:)` / `idx(:)`.
- Dead code: commented `snapclim_end`, the inline commented temperature-correction block, the
  empty "write static fields" stub, the commented type-dump TODO.

## Status / next steps

- **Done:** design agreed; worktree `snapclim2` created; compiling module skeleton
  (`libs/snapclim2.f90`) with the type layout and stub public API; Makefile rule.
- **Next:** implement `snapclim2_par_load` (top-level + registry + snapshot + index groups),
  then the load→combine→transform pipeline, then the validation harness against snapclim.

## Open items

- Final module/type name.
- Exact registry field set to name explicitly vs leave to `extra(:)`.
- Whether `dTa/dTo/dSo` (driver tsforcing) and snapclim2's own `idx(:)` should compose or be
  mutually exclusive per field (today snapclim uses the arg when present, else its own index).
- Provenance record format (reuse ncio; decide whether derived-field snapshot is default-on).
