# CISM as an ice-sheet backend: evaluation and plan

Status: evaluation only (2026-10-04). The leading `_` keeps
quarto from rendering it. To be continued in a dedicated new repo for the
shared ice-model interface.

States checked: yelmox dev e117800, climber-x (local checkout), CISM 71291a4
(tag `cism_main_3.00.001`, ISMIP7 merged).

## Verdict

Feasible. The right route is a shared `ice_model` interface library used by
both yelmox and CLIMBER-X, with Yelmo and CISM as backends (SICOPOLIS stays
in CLIMBER-X). Most of the work is in yelmox, which has no layer between the
driver and Yelmo. CISM can be embedded with a thin wrapper around glissade;
none of its existing coupling APIs fit as they are.

## yelmox (dev)

- No isolation layer: `kryos_domain` holds `type(yelmo_class) :: yelmo`
  (`libs/kryos.f90:144`). Nearly every kryos object is built with
  `$(INC_YELMO)`, and many modules take `wp` from `yelmo`
  (e.g. `kryos_forcing.f90:12`).
- Seams that already exist:
  - `couple_*_to_yelmo` (`libs/kryos_coupling.f90:131-315`): all boundary
    condition writes.
  - `couple_yelmo_to_htopo` (`:433`): geometry reads for the hub.
  - `step_icesheet` (`:317`), which calls `yelmo_update` (`:330`).
  - Isostasy, marine shelf, SMB and climate only see the hub (`dom%topo`),
    never Yelmo.
- Boundary conditions written into Yelmo:
  - isostasy and sea level: `bnd%z_bed = z_bed_ref + remap(w+we)`,
    `bnd%z_sl = bsl + remap(z_ss-bsl)` (coupling:196-197)
  - SMB: `bnd%smb` (via `conv_we_ie*1e-3`) and `bnd%T_srf` (:279-280), plus
    optional glacial scaling and a `lim_pd_ice` mask (:284-295)
  - marine shelf: `bnd%bmb_shlf`, `bnd%T_shlf` (:313-314)
  - climate: `bnd%Qd` (:154); geothermal: `bnd%Q_geo` (kryos.f90:288);
    sediments: `bnd%H_sed` (:284)
  - basal hydrology is internal to Yelmo; the driver writes no hydrology
    fields
- Fields read back from Yelmo:
  - hub: `tpo%now%H_ice, f_grnd, z_srf`, `bnd%z_bed, z_bed_ref, z_sl`
    (coupling:450-459)
  - isostasy: `H_ice`, `bnd%dzbdt_corr` (:117-118); at init `bnd%z_bed_ref`,
    `H_ice_ref` (startup:200-204)
  - forcing feedback: `reg%V_ice` (yelmox.f90:156)
  - bipolar FWF: `tpo%now%mb_net/smb/bmb/cmb/H_ice/dHidt/f_grnd`
    (obm_coupling.f90:199-205)
  - CMIP output: about 30 `tpo`/`dyn`/`thrm` fields (cmip_output.f90)
- Code tied to Yelmo internals (the hard part of a refactor):
  - startup and restart: `yelmo_init_state`, `yelmo_update_equil`,
    `yelmo_init_topo`, `yelmo_restart_read_topo_bnd`, `yelmo_restart_init`,
    `calc_ytopo_diagnostic`, `yelmo_regions_update` (kryos_startup)
  - cold-start edits to `tpo%now%H_ice`, `bnd%H_ice_ref`, `bnd%mask_ice`
  - opt tuning: `tpo%par%topo_rel*`, `dyn%now%cb_ref`, `dyn%par%till_*`,
    `calc_cb_ref`, `ice_optimization`
  - NEGIS `cb_ref` overwrite (kryos_regions:104-116)
  - regions, ts output and CMIP output (`yelmo_regions_*`,
    `yelmo_write_reg_init`)
- Restart bundle (`restart-<kyr>-kyr/`): `yelmo_restart.nc` sits alongside
  the isostasy, shelf, smbpal, climate, bsl and tsforcing files. The hub is
  rebuilt from the restored models.
- Pattern to copy: the climate backend (`libs/yelmox_climate.f90`). It uses
  a runtime string (`[coupling] climate`), `select case` dispatch, a neutral
  output type (`climate_out_class`) with capability flags, and an optional
  library behind a stub (`climate_rembo.f90` / `climate_rembo_stub.f90`,
  `rembo=1`).

## CLIMBER-X `src/ice/`

- `ice_def.f90:39-76`: a backend-neutral `ice_class` with SI units, fields on
  the ice grid, and a lookup from the ice grid to the coupler grid.
  - Inputs: `mask_ocn_lake, z_sl, z_bed, z_bed_fil, smb, accum, runoff,
    bmb_shlf, temp_s, temp_g, q_geo, H_sed, t_ocn, s_ocn`.
  - Outputs: `mask_extent, z_sur, z_sur_std, z_bed_std, z_base, H_ice, calv,
    Q_b, error`.
- `ice_model.f90`:
  - Public procedures: `ice_init_domains, ice_init, ice_update, ice_end,
    ice_write_restart`.
  - Dispatch is a runtime `select case(model)` over the module-global
    `ylmo_doms(:)` and `sico_doms(:)` (`:46-47`).
  - Both libraries are always `use`d (`:35-39`), and there is no stub per
    backend. `.ice_model_dummy.f90` replaces the whole module in climate-only
    builds.
- Yelmo wrapper: `ice_to_yelmo` (`:534-569`) and `yelmo_to_ice`
  (`:494-532`). `z_sl = 0` is hard-coded, and accum, runoff, t_ocn, s_ocn and
  mask_ocn_lake are unused.
- SICOPOLIS: an in-tree port (`src/ice_sico/`) with (j,i), 0-based arrays.
  The wrapper is at `ice_model.f90:571-729`.
- Gaps:
  - no velocity or mass-flux outputs
  - restarts differ by backend (Yelmo uses `restart=` in its own namelist and
    ignores `l_restart`)
  - output is not part of the interface
  - `ice_end` is never called
- Regridding happens in `src/main/coupler.f90`:
  - geo to ice: conservative mapping, nearest-neighbour for masks
  - ice to coupler grid: calv and Q_b summed per cell
  - smb and bmb already run on the ice grid
- Multiple domains: up to 10 named domains, looped serially. The ice model is
  called every `n_year_ice` years.

## CISM

### Build

- CMake builds one static library, `libglimmercismfortran.a`
  (CMakeLists.txt:433), plus the generated `.mod` files in
  `${CISM_BINARY_DIR}/include`. `CISM_BUILD_CISM_DRIVER=OFF` builds only the
  library.
- Settings for embedding: `CISM_SERIAL_MODE` (no MPI; `parallel_slap.F90`),
  Trilinos off, `CISM_USE_GPTL_INSTRUMENTATION=OFF`.
- NetCDF is required (`CISM_NETCDF_DIR`). Python is needed at build time
  (`utils/build/generate_ncvars.py` generates `glide_io.F90` etc.).

### Existing APIs and why they don't fit

- `cism_driver/` (`cism_front_end.F90`): owns the whole time loop and reads
  the config path from the command line.
- GLINT: built for a GCM on a global lat-lon grid, with downscaling and PDD.
- GLAD (`libglad/glad_main.F90`): fields on the CISM grid, but inputs are
  limited to SMB, T_srf and ocean T/S. Time is in integer hours with a fixed
  accumulate/average cadence. topg, sea level, Q_geo, bmb and thickness
  cannot be passed in.

### Proposed wrapper (modelled on `glad_initialise.F90` and `glad_timestep.F90`)

- Init:
  - `ConfigRead`, then override values with `ConfigSetValue`
  - `glide_config`
  - `glissade_initialise(model, evolve_ice)` (glissade.F90:89)
  - `glissade_diagnostic_variable_solve` (:2109)
- Step: set the `model%...` fields, then call `glissade_tstep(model, time)`
  with time in years (:1172).
- End: `glide_finalise` (glide_stop.F90).

### Where fields live (`libglide/glide_types.F90`)

| Field | Location | Notes |
|---|---|---|
| `thck`, `usrf`, `lsrf`, `topg`, `f_ground`, masks | `model%geometry` (1201-1262) | |
| `acab` | `model%climate%acab` (1526) | m ice/s |
| `smb` | `model%climate%smb` (1531) | mm/yr w.e., used via `smb_input` |
| `artm` | `model%climate%artm` (1534) | |
| `eus` | `model%climate%eus` (1586) | scalar sea level |
| `bheatflx` | `model%temper%bheatflx` (1788) | |
| `bmlt_float_external` | `model%basal_melt` (1898-1907) | used with `whichbmlt_float=4` (external field) |
| `thermal_forcing(nzocn,:,:)` | `model%ocean_data` | used with `whichbmlt_float=6` (melt from thermal forcing) |
| calving and mass fluxes | `model%calving`, `model%mass_flux` | |
| velocities | `model%velocity` | |

### External bed and sea level

- With `isostasy = 0`, CISM never changes `topg`.
- Caveat: the topg halo update and the `marine_connection_mask` recompute
  only run inside the ISOSTASY_COMPUTE branch
  (glissade_isostasy.F90:273-310; see comment at glissade.F90:896). The
  wrapper must do both after setting topg or eus. This may need a small CISM
  patch or calls to those internals.

### Grid and parallelism

- Regular x-y grid (`[grid]` ewn, nsn, dew, dns), decomposed into blocks with
  a 2-cell halo.
- Serial mode: no MPI, communicator variables are constants.
- With MPI, a host model would call `parallel_set_info(comm, main_rank)`
  instead of `parallel_initialise`, and use `gather_var` / `scatter_var` and
  the haloed/non-haloed conversions.

### Configuration and I/O

- Config is INI-style `.config` files, read only from a file
  (`glimmer_config.F90:105`).
- CISM always does its own NetCDF input (`glide_io_readall`), output
  (`glide_io_writeall`) and `[CF restart]` files (`restart` = 0, 1 or 2).

### Obstacles

- A fatal error calls `parallel_stop`, i.e. `mpi_abort` and `stop`; there are
  also about 20 bare `stop` statements. The driver cannot recover from these.
- Module-global state: the communicator, `nhalo`, and a model registry
  (max 8). Multiple instances share one communicator.
- Some options are not parallel-safe (e.g. GTHF_COMPUTE).
- Performance: yelmox is OpenMP-only. Serial CISM may be slow at high
  resolution, and MPI CISM inside yelmox needs extra wrapper work.

## Plan (separate steps, each verified on Levante)

1. **Shared interface library (new repo)**
   - Neutral `ice_class` covering what yelmox and CLIMBER-X exchange:
     - inputs: z_bed, z_sl (2D), smb, T_srf, bmb_shlf, T_shlf / t_ocn /
       s_ocn, Q_geo, H_sed, Qd, mask_ice, regions and basins
     - outputs: H_ice, z_srf, z_base, f_grnd, masks, dHidt,
       smb / bmb / cmb / mb_net, calving, ux_s / uy_s, z_bed_ref /
       H_ice_ref, dzbdt_corr
   - Procedures: `ice_init, ice_update(time), ice_restart_write/read,
     ice_write_* (handed to the backend), ice_end`.
   - Yelmo backend: merge CLIMBER-X `ice_to_yelmo` / `yelmo_to_ice` with
     yelmox `couple_*_to_yelmo` and `couple_yelmo_to_htopo`.
2. **yelmox on the interface, Yelmo only**
   - Replace `dom%yelmo` and keep the restart bundle layout.
   - Goal: bit-identical to the current baseline. Worth doing even without
     CISM.
3. **CISM backend**
   - Build flag `cism=0|1` with a stub module (rembo pattern).
   - Generate the CISM `[grid]` and options from `[domain]` via
     `ConfigSetValue`, so there is one source of truth for the grid.
   - Convert units (acab in m ice/s, T in K vs degC).
   - topg/eus halo update and mask recompute after each external change.
   - The CISM restart file becomes part of the restart bundle.
   - First tests: serial, ANT-32KM and GRL-16KM, compared with Yelmo.
4. **CLIMBER-X**
   - Replace `src/ice/ice_model.f90` with the library.
   - SICOPOLIS stays a backend defined in CLIMBER-X.
   - Check that results are identical with Yelmo.
5. **Later:** CISM with MPI, and tuning CISM physics for the experiments
   (likely the biggest scientific effort).

## Open decisions

1. **Where the interface lives:** a new standalone repo (intended), inside
   yelmo, or inside fesm-utils?
2. **How backends are selected:**
   - An abstract type that each backend extends. A host program can add its
     own backend without the library knowing it, which CLIMBER-X needs for
     SICOPOLIS. Recommended.
   - Or `select case` plus stubs, matching the yelmox climate backend.
3. **Yelmo-only features** (opt / `cb_ref`, NEGIS, equilibrium spin-up, CMIP
   and regions output):
   - restrict them to the Yelmo backend through a typed accessor
   - or generalise them as capabilities (e.g. CISM's own basal inversion)
4. **CISM in serial first,** leaving MPI for later?
5. **Output:** each backend writes its own files (CISM does anyway), or the
   driver also writes a neutral ice file for analysis that works the same
   for every backend?

## Addendum: how CISM is coupled in CESM (2026-10-04)

Sources: ESCOMP/CISM-wrapper main@718f8f27 (the GLC cap), ESCOMP/CMEPS
main@869eed8c (the mediator), and local CISM `libglad/`.

### Parallelism

- Grids stay distributed end to end; nothing is gathered on the coupling path.
- GLC runs on its own PEs. The cap passes the component communicator to CISM
  via `parallel_set_info` (`mpi/glc_communicate.F90`).
- `glc_indexing.F90` builds the ESMF mesh and DistGrid from CISM's owned
  (non-halo) points:
  - local sizes from `glad_get_grid_size`, positions from
    `glad_get_grid_indices`
  - global index = `(row-1)*global_ewn + col`
- Fields go MED -> GLC by redist. The mediator regrids with ESMF route
  handles: lnd->glc and ocn->glc bilinear, glc->lnd and glc->rof
  conservative.
- libglad only touches local points: haloed <-> non-haloed copies plus a
  neighbour halo exchange. Gathers happen only in CISM's own netCDF I/O.
- Serial `parallel_slap` has the same API.

### Time coupling

- CISM is called daily (`GLC_NCPL=1`/day). The mediator accumulates lnd and
  ocn inputs and averages them yearly (`GLC_AVG_PERIOD`).
- glad's mass-balance step is hard-coded to 1 yr (`glad_mbal_coupling.F90:101`).
- `valid_inputs` must be true at the end of that step, or glad aborts.
- The ice dycore then runs `ice_tstep_multiply * mbal_accum_time / ice_tstep`
  steps.
- Time is integer hours on a 365-day calendar. The driver owns the clock.

### Fields and processing

- **Imports:** qsmb (kg m-2 s-1 w.e.), tsfc (degC), ocean T (K) / S on 30
  levels.
- **Ocean coupling is not actually wired:** the cap never unpacks T/S and
  hard-codes 274 K / 35 (`glc_InitMod.F90:398-404`, TODO).
  `ocean_data_domain=2` through glad is effectively untested.
- **SMB downscaling happens in the mediator, not CISM:**
  - CTSM sends SMB in elevation classes.
  - These are interpolated in elevation using CISM `usrf`.
  - Optional global renormalisation (`glc_renormalize_smb`) uses separate
    factors for accumulation and ablation.
- **Areas:** `glad_get_areas` returns the projection area dx*dy with no map
  factor. CMEPS deliberately applies no area correction for glc.
- **Exports:**
  - `icemask` (usrf > 0), `ice_covered` (thck > 0), `topo` (usrf in the
    mask)
  - `rofi` (calving + removal), `rofl` (basal melt), `hflx`
  - Fluxes are averages over the last mass-balance period and are zeroed
    unless `zero_gcm_fluxes=.false.` (two-way and evolving).
- **No bed, isostasy or sea-level exchange.** Land ice topography reaches the
  atmosphere only through CTSM elevation classes.

### Multiple sheets and restarts

- One GLC component holds N glad instances (GrIS, AIS), each with its own
  config and mesh. All share the same mass-balance and ice timesteps
  (`check_mbts`).
- The cap writes one restart per sheet. Partial glad input averages are not
  saved: `glad_okay_to_restart` requires `av_steps == 0`. Mid-year restarts
  work only because the mediator holds the accumulators.

### Relevance to our effort

- **Parallelism (no issue):** we run single-process, so CISM can be serial or
  MPI with one rank (`parallel_set_info(MPI_COMM_SELF, 0)`). Then own =
  global, and glad arrays are plain `(nx, ny)` with halos hidden. The CESM
  pattern (mesh from owned points + global indices) shows how to go parallel
  later.
- **Does glad suffice? No:** it has no inputs for topg, z_sl/eus, bmb_shlf or
  Q_geo, and CESM never exercises ocean forcing through it. The topg/eus halo
  and marine-mask issue is ours alone, because CESM never changes the bed
  from outside.
- **Our wrapper:** model it on glad but extend it. Two options:
  - (a) our own thin glissade wrapper, as planned
  - (b) extend glad upstream with optional topg/eus/bmb/Q_geo inputs and
    contribute that back to ESCOMP/CISM
- **SMB is our job:** downscaling and conservation are the caller's
  responsibility. yelmox already computes SMB on the ice grid at the current
  surface, and CLIMBER-X does too. Renormalisation would be the driver's
  choice.
- **Areas:** CISM integrates in projection area dx*dy. Any integral in true
  area (sea-level contribution, freshwater to CLIMBER-X `ice_to_cmn`) needs a
  map-factor correction done by us. The convention must match Yelmo's so
  backends are comparable.
- **Flux semantics map directly:**
  - rofi (calving + removal) -> CLIMBER-X `calv`
  - rofl (basal melt) -> `Q_b`
  - glad `ice_covered` / `icemask` correspond to Yelmo `mask_ice` /
    `H_ice > 0`
- **Restarts:** take them only at the end of a mass-balance period. That is
  naturally satisfied if the ice backend is called every 1 yr (or more).
- **Multiple domains:** N instances in one process matches yelmox_bipolar and
  CLIMBER-X `ice_domain_name(:)`. CISM limits are a registry max of 8 models
  and a single shared communicator.
- **Not verified:** CISM's own OpenMP coverage; how CTSM consumes glc fields.

### Additional open decision

6. **CISM wrapper:** our own glissade wrapper, or extend glad upstream?
