---
title: "yelmox"
aliases:
  - flavor-yelmox.html
  - flavors.html
---

The YelmoX program: one ice-sheet domain with FastIsostasy bedrock, a shared
barystatic sea level, a climate backend chosen at runtime, a surface mass balance
and the marine shelf, built on the multigrid `kryos_domain` (see
[Multigrid coupling](multigrid.md)). Its variant
[`yelmox_bipolar`](yelmox-bipolar.md) runs two domains, north and south, coupled
through the sea level and an ocean box model.

| Program | Build | Domains | Climate | Ocean |
|---|---|---|---|---|
| `yelmox` | `make yelmox` (`rembo=1` for REMBO) | one | any backend (`[comps] climate`) | the climate (depth profiles or shelf base) |
| [`yelmox_bipolar`](yelmox-bipolar.md) | `make yelmox_bipolar` | north + south | snapshots (×2) | snapshots + shared ocean box model |

- **Program:** `yelmox/yelmox.f90` (thin driver) + `libs/kryos*.f90` (domain, coupling, startup, output).
- **Build:** `make yelmox` (`make yelmox rembo=1` to link REMBO).
- **Configs:** `yelmox/yelmox_<domain>.nml` (Antarctica, Greenland, North, LIS,
  Pyrenees, SRG, plus `pd_` present-day and paleo variants), `yelmox_esm_*.nml`
  ([ESM forcing](climate-esm.md)), `yelmox_rembo_Greenland.nml`
  ([REMBO](climate-rembo.md)), `yelmox_Greenland_ramp.nml` (an air temperature
  ramp: `[anom]` on a `ramp-slope` index).

## Components

| Role | Module | Grid |
|---|---|---|
| Ice sheet | Yelmo | `grid_ice` |
| Isostasy + sea level | FastIsostasy (`isos`) + shared `bsl` | `grid_isos` |
| Climate (atmosphere + ocean) | backend of `[comps] climate`: [snapshots](climate-snapshots.md), [esm](climate-esm.md) or [rembo](climate-rembo.md) | `grid_clim` |
| Surface mass balance | chion (default), smbpal, `smb_simple` or the climate's own | `grid_surface` |
| Sub-shelf melt | marine_shelf | `grid_shelf` |
| Geometry hub | htopo | `grid_hub` (hi-res) |

Each module runs on its own grid, set in `[domain]`; the coupler remaps fields between
grids at the moment of coupling. See [Multigrid coupling](multigrid.md).

### Climate backends

The forcing of a run is set at runtime by its climate backend, `[comps] climate`:

| Backend | Supplies | Page |
|---|---|---|
| `snapshots` | atmosphere + ocean profiles, from climate snapshots blended by indices, one blend model over any number of snapshots | [Snapshots](climate-snapshots.md) |
| `esm` | reference climatology + Earth-system-model anomalies; atmosphere or smb, ocean at the shelf base, subglacial discharge | [ESM forcing](climate-esm.md) |
| `rembo` | REMBOv1 atmosphere + smb; ocean from snapshots (ocean-only snapshot) | [REMBO climate](climate-rembo.md) |

## Stepping order

The driver owns the timeline (`ts`) and the shared sea level (`bsl`), and advances
the domain once per step with the coupling sequence written out in the time loop:

```fortran
call step_relax(dom, ts)          ! topography relaxation (relax)
call step_optimize(dom, ts)       ! cb_ref/tf_corr optimization (opt)
call step_isostasy(dom, ts, bsl)  ! bedrock + sea level (dt_isos cadence)
call couple_to_yelmo(dom)         ! bedrock now; smb + shelf melt lag one step
call step_icesheet(dom, ts)       ! yelmo_update
call couple_yelmo_to_htopo(dom)   ! hi-res geometry from the models
call step_climate(dom, ts, tsf)   ! climate (dt_clim cadence)
call step_surface(dom, ts)        ! surface mass balance + temperature
call step_shelf(dom, ts)          ! shelf-base melt + temperature
```

`couple_to_yelmo` assembles the Yelmo boundary state: isostasy from this step,
smb and shelf melt from the **previous** step (a one-step coupling lag);
`step_climate`, `step_surface` and `step_shelf` then produce the forcing
consumed on the next step.
The climate is refreshed on the `comps.dt_clim` cadence, isostasy on the
`comps.dt_isos` cadence (the bedrock is held in between); the smb every step.

## Indices and the homogeneous anomaly (`[indices]`, `[anom]`)

**Indices** are named scalar series, each a `tsgen` time-series generator (the
modern replacement for the legacy `hyster` module). They are **driver-owned**:
the program holds one registry (`indices_class`, `idx`, libs/kryos_forcing.f90),
advances every index once per step on the forcing clock `time_rel`, and passes
it to `step_climate` (and the cold start); consumers read an index by name
(the snapshots axes, `[anom]`). Multi-domain drivers share one registry.

```fortran
&indices
    names = "ramp"            ! "" = no indices
/

&idx_ramp                     ! a tsgen group (unset keys: template [idx])
    method    = "ramp-time"   ! const | ramp-slope | ramp-time | ramp-time-step | sin | series | exp | PI42 | ...
    f_min     = 0.0
    f_max     = 5.0
    dt_ramp   = 200.0
    ! ... (see the [idx] template of input/yelmox_defaults.nml)
    ! units = "K" (label of idx_<name> in the 1D output, default "1"),
    ! restart_every_df = 0.5 (a restart bundle each |Δf|, default 0 = off)
/

&anom
    active = True     ! a homogeneous anomaly from one index
    index  = "ramp"
    f_ta   = 1.0      ! dTa = f * f_ta   (air temperature [K])
    f_to   = 0.0      ! dTo = f * f_to   (ocean temperature [K])
    f_so   = 0.0      ! dSo = f * f_so   (ocean salinity [psu])
/
```

Time-driven methods (`ramp-*`, `sin`, `const`, `series`) are analytic or
tabulated; feedback methods (`exp`, PI/PID controllers) modulate the forcing
rate from the model response — the response variable passed to `tsgen` is
total ice volume (Gt). The run stops when any index trips its kill switch
(`with_kill`), and writes a restart bundle `restart-<n>` each time an index
moves by more than its `restart_every_df`. Each index restarts from
`idx_<name>_restart.nc` in the bundle, and is written to the 1D output as
`idx_<name>` (plus `idx_<name>_df_dt`).

`[anom]` is a climate modifier (libs/climate_mods.f90): applied to the
backend's output after the backend, before `[ocn_fraction]`, `[ocn_obm]` and
`[bias]`, on the climate grid. Air temperature + `dTa`, precipitation × exp(β·`dTa`) (β the backend's
precipitation sensitivity), ocean profiles + `dTo`/`dSo` (esm: at the shelf
base). The group is optional (inactive without it) and takes the domain suffix
(`anom<sfx>`). Checked at init: the index must be in the registry, and a
backend that supplies the smb (esm with `surface_method = "climate"`) refuses
it. rembo reads `[anom]` itself, as an input of REMBO (see
[REMBO](climate-rembo.md)).

## Ocean anomaly from the atmosphere (`[ocn_fraction]`)

The climate modifier after `[anom]`: the ocean temperature profiles, at every
depth, are the reference plus `f_to` times the mean sea-level air temperature
anomaly of the current update over a region (salinity unchanged):

    to_ann = to_ann_ref + f_to * mean_region(tsl_ann - tsl_ann_ref)

The anomaly is that of the same update (no lag) and includes `[anom]`'s
`dTa`; at sea level, it excludes the lapse-rate cooling of the ice sheet.

```
&ocn_fraction
    active = True
    f_to   = 0.25     ! [1] fraction of the air temperature anomaly
    region = "all"    ! cells of the mean
/
```

`region` is `"all"` or a selection expression of the domain's regions on the
climate grid (`[regions]`, resolved once at init), e.g. `"zone:open_ocean"`,
`"zone:open_ocean,continental_shelf"`, `"region:Greenland"`, `"~zone:land"`
(terms `layer:values` joined by `&` and `|`; see the fesm-utils regions docs).
Its cell count is printed at init; an empty region stops the run. Refused for
a backend without ocean profiles or sea-level air temperature (esm, rembo),
and with `[anom] f_to /= 0` (both would set the ocean temperature anomaly).
The anomaly is written to the climate's 2D file (`ocn_fraction_dto`). See
[Snapshots](climate-snapshots.md).

## Ocean from the ocean box model (`[ocn_obm]`)

The climate modifier after `[ocn_fraction]`, for
[`yelmox_bipolar`](yelmox-bipolar.md): the ocean temperature profiles, at every
depth, from the box of the domain's hemisphere of the ocean box model
(northern box `tn`, southern box `ts`); salinity unchanged:

    absolute:  to_ann = box
    anomaly:   to_ann = to_ann_ref + (box - box_init)

```
&ocn_obm_north
    active = True
    mode   = "absolute"   ! absolute | anomaly
/
```

The one modifier not applied in the climate step: the driver that owns the box
model applies it (`climate_mod_ocn_obm`) after both domains' climate steps and
before their shelf steps; `active` is the switch of that exchange
(`yelmox_bipolar` needs `[ctrl] active_obm`). `box_init` is the box's initial
value of the par file (`tn_init`, `ts_init`), also after a restart. Refused for
an ocean at the shelf base (esm), and with `[ocn_fraction]` or `[anom] f_to /= 0`
(it replaces the ocean temperature). The box temperature and its anomaly are
written to the climate's 2D file (`ocn_obm_theta`, `ocn_obm_dto`).

## Climate bias correction (`[bias]`)

A static spatial correction of the atmosphere, the last climate modifier,
applied on the climate grid to the monthly air temperature and precipitation
of any backend except rembo (annual fields only), before the move to the
surface grid's elevation:

    tas = tas + f_tas*dtas,    pr = pr*(1 + f_pr*(fpr - 1))

The annual means follow (temperature additively, precipitation by the ratio
of the annual totals); the reference climate and the ocean are unchanged. The
group is optional; without it, or with `active = False`, nothing changes.

```
&bias
    active   = True
    var_defs = ""             ! files with the groups below, searched in order, then this file
    dtas     = "bias_dtas"    ! var_defs group: temperature correction [K]; "" = none
    fpr      = "bias_fpr"     ! var_defs group: precipitation factor [1]; "" = none
    f_tas    = 1.0            ! weight of dtas
    f_pr     = 1.0            ! weight of fpr
/
```

`dtas` and `fpr` are varslice groups (as snapshots `var_defs`): one field for all
months (`time = 0 0 0 0 1`) or 12 monthly fields (`time = 0 0 0 0 12`), on the
climate grid or remapped online (`remap`, `grid_src`). The weighted annual
fields are written to the climate's 2D file (`bias_dtas_ann`, `bias_fpr_ann`).

## Insolation and orbit (`[orbit]`)

The surface models that take the top-of-atmosphere insolation (chion's ITM and
BESSI, smbpal's ITM) compute it with the [insol](https://github.com/fesmc/insol)
library (a dependency, `libinsol.a`), the daily mean in closed form from the
orbital elements of Laskar et al. (2004). `[orbit]` (one per domain,
`[orbit_north]` in yelmox_bipolar) gives the elements; both surface models read
it:

```
&orbit
    method = "time"   ! clock: at the forcing clock time_rel; time: at a fixed time; fixed: ecc, obl, per
    time   = 0.0      ! [yr] model time of the orbit (time; relative to 1950, past < 0)
    ecc    = 0.016702 ! [1]   eccentricity (fixed)
    obl    = 23.4393  ! [deg] obliquity (fixed)
    per    = 102.918  ! [deg] longitude of perihelion (fixed)
    fldr   = "input"  ! folder of the orbital tables
/
```

`clock` follows the forcing clock (`ts%time_rel`: transient insolation with
`tstep_method = "rel"`); `time` holds the orbit of one time (present day by
default, the 1950 orbit); `fixed` takes the given elements (the defaults are the
tables' at 1950). A fixed orbit is the one exception to the single forcing
clock. Each model recomputes its insolation only when the elements change, so
once per run unless `method = "clock"`. The par files list `method` and `time`;
the other keys come from `input/yelmox_defaults.nml`.

## Configuration

A run is set by these namelist groups. Every key is required (a missing key stops
the run), except the keys of an unselected `init_method` and the `[relax]` and
`[opt]` groups when their switch in `[sim]` is off.

The forcing groups have defaults: `input/yelmox_defaults.nml` declares every
key of `[indices]`, `[idx_<name>]` (template `[idx]`: the tsgen keys, `units`,
`restart_every_df`), the snapshot backend (`[snap]`, templates `[snap_field]`,
`[snap_snapshot]`, `[snap_axis]`), the climate modifiers (`[anom]`,
`[ocn_fraction]`, `[ocn_obm]`, `[bias]`), the orbit of the insolation (`[orbit]`)
and `[obm_forcing]` (yelmox_bipolar),
with the domain suffix of a multi-domain driver
(`[anom_north]` to `[anom]`). A par file lists only what differs; an absent
group takes every default (a modifier is then inactive). Each read goes through
`nml_read(..., defaults_file, defaults_group)` (the path in `libs/yelmox_nml.f90`),
and a key not declared there stops the run (`nml_validate`).
`scripts/check_par_nml.py` checks the yelmo groups against
`input/yelmo_defaults.nml` and these groups against `input/yelmox_defaults.nml`
offline (`--dead`: keys of the other groups that no source reads). The other
yelmox groups below are read without defaults.

| Group | Sets |
|---|---|
| `[ctrl]` | the run: its timeline group, calendar years, the restart bundle to start from |
| timeline group (`[ctrl] run_step`) | the timeline: `tstep_method`, `tstep_const`, `time_init`, `time_end`, `dtt` |
| `[domain]` | the domain: name, grid of every component, hub topography and masks (see [Multigrid coupling](multigrid.md)) |
| `[comps]` | the components: which are active, with which model, how often |
| `[sim]` | the conditions of this simulation: cold-start ice state, relaxation, optimization, regional modifications |
| `[relax]`, `[opt]` | the parameters of the relaxation and the optimization switched on in `[sim]` |
| `[output]` | the output files; their intervals are `[tm_1D]`, `[tm_2D]`, `[tm_2Dsm]`, and `[tm_rst]` for restarts |

The models have their own groups (`[yelmo]` and the Yelmo physics groups,
`[isos]`, `[barysealevel]`, `[marine_shelf]`, `[surface_chion]` + `[chion]`, `[smbpal]`, the climate backend's,
...). In `yelmox_bipolar` the domain groups carry the hemisphere suffix
(`[domain_north]`, `[comps_north]`, `[sim_north]`, `[relax_north]`,
`[opt_north]`, `[output_north]`, ...); `[ctrl]` and the timeline are shared.

A minimal skeleton:

```fortran
&ctrl
    run_step        = "ctrl"            ! timeline group ("ctrl" = this group)
    calendar        = False
    calendar_ref    = 1950.0
    restart         = "None"            ! restart bundle folder, or "None" for a cold start
    tstep_method    = "const"
    tstep_const     = 0.0
    time_init       = 0.0
    time_end        = 15e3
    dtt             = 10.0
/

&comps
    with_ice_sheet  = True
    with_isostasy   = True
    with_climate    = True
    with_surface    = True
    with_shelf      = True
    climate         = "snapshots"
    surface_method  = "chion"
    dt_clim         = 10.0
    dt_isos         = 10.0
/

&sim
    init_method       = "none"
    init_marine_H     = False
    init_kill_shelves = False
    init_time_thrm    = 0.0
    relax             = True            ! uses [relax]
    opt               = True            ! uses [opt]
    use_negis         = False
    scale_glacial_smb = False
    lim_pd_ice        = False
/
```

### `[ctrl]`

| Key | Values | |
|---|---|---|
| `run_step` | group name | the group holding the timeline: `"ctrl"` for `[ctrl]` itself, or a phase group such as `[spinup]` or `[transient]`, so one par file can hold several run phases |
| `calendar`, `calendar_ref` | bool, [yr CE] | timeline in calendar years, against the reference year `calendar_ref` (ESM runs) |
| `restart` | folder, `"None"` | the restart bundle to start from; `"None"` = cold start |

The timeline group holds `tstep_method` (`"const"`: the time advances and the
forcing time stays at `tstep_const`; `"rel"`: the forcing follows the time, e.g.
for paleo records; `"cal"`: calendar time), `tstep_const`, `time_init`,
`time_end` and the main time step `dtt` [yr].

A restart bundle is a folder `restart-<kyr>-kyr/`, written on the `[tm_rst]`
schedule and always at `time_end`. It holds one restart file per stateful model,
the shared sea level and the transient-forcing state (in `yelmox_bipolar` also the
ocean box model, with each domain in a subfolder named after it). Pass `restart`
as an absolute path: the run starts inside its own folder.

### `[comps]`

The components: which are active, with which model, how often.

| Key | Values | |
|---|---|---|
| `with_ice_sheet`, `with_isostasy`, `with_climate`, `with_surface`, `with_shelf` | bool | components in the coupling sequence: the ice sheet, isostasy, the climate (atmosphere + ocean), the surface (mass balance + temperature) and the shelf base (melt + temperature); `with_surface` and `with_shelf` need `with_climate` |
| `climate` | `snapshots`, `esm`, `rembo` | the climate backend |
| `surface_method` | `chion` (default), `smbpal`, `smb_simple`, `climate` | chion or smbpal (ITM or PDD), both from the climate's temperature and precipitation; chion runs the ITM, BESSI or PDD snowpack (`[chion] model`) in daily steps (`[surface_chion]`, `[chion]`; ITM's parameters are the group `[chion] nml_itm` names: smbpal's `[itm]`, or for Greenland `[itm_chion]` (`alb_ice = 0.31`, calibrated to MAR); BESSI gets the surface shortwave as TOA insolation times `trans_sw`, the TOA insolation itself (chion's cloud-proxy longwave), the ice thickness, a constant `wind_speed` and `rel_hum`, and the air pressure from the surface elevation; its constants are `[chion_const]`, for Greenland `alpha_ice = 0.55`, calibrated to MAR; BESSI runs pass `surface_chion.time_equil=20`, converged within 0.3 Gt/yr of the par files' 100 yr, which ITM needs; PDD's parameters are the group `[chion] nml_pdd` names, `[pdd_chion]` with `pdd_method = "cg05"` (Calov and Greve, 2005), smbpal's PDD; the surface temperature of ITM and PDD is chion's `chion_tsrf_annual`, the annual-mean air temperature warmed by `firn_fac` times the annual refreezing, capped at freezing over ice); smb_simple (needs a sea-level air temperature: snapshots); the climate's own smb (esm, rembo; required by rembo) |
| `dt_clim` | [yr] | climate update interval; `<= 0`: updated only at the cold start |
| `dt_isos` | [yr] | isostasy update interval (> 0; 10 in all par files): FastIsostasy integrates the whole interval in one call, the bedrock is held in between |

### `[sim]`

The conditions of this simulation: the cold-start ice state (`init_*`), the
relaxation, the optimization and regional modifications. The `init_*` keys act on
a cold start only; `relax` and `opt` act on a cold start or a restart, with their
times counted from the start of the run (`time_init`).

| Key | Values | |
|---|---|---|
| `init_method` | `none`, `equil`, `recon`, `recon_ref` | cold-start ice state: as initialized; a short equilibration with constant boundaries; the reconstruction `recon_path` as initial ice where `recon_regions` selects; the reconstruction as reference ice only |
| `init_equil_time` | [yr] | `equil`: equilibration time |
| `recon_path`, `recon_var` | path, name | `recon`, `recon_ref`: the reconstruction file (`{domain}`, `{grid_name}` = `grid_ice`) and its ice-thickness variable |
| `recon_time` | yr | `recon`, `recon_ref`: model time of the state; the nearest slice of `recon_path` (its time axis converted from the units: `years ago`/`BP` count into the past, `kiloyears`/`ka` in thousands) |
| `recon_regions` | selection expression | `recon`: where its ice is imposed (e.g. `"region:North_America"`) |
| `init_marine_H` | bool | LGM-like marine ice, before `init_method` |
| `init_kill_shelves` | bool | no ice where the present-day bed is ocean |
| `init_time_thrm` | [yr] | then equilibrate with the topography fixed (`0` = off) |
| `relax` | bool | relax the topography towards the reference (`[relax]`) |
| `opt` | bool | optimize the basal friction and the thermal-forcing correction (`[opt]`) |
| `use_negis` | bool | NEGIS basal-friction modification (`[negis]`) |
| `scale_glacial_smb` | bool | reduce negative glacial smb (`[glacial_smb]`) |
| `lim_pd_ice` | bool | extra melt (4 m/yr) outside the present-day ice extent |

A spin-up usually sets `relax` and `opt` together: the relaxation holds the
floating ice and the grounding zone near the observations while the optimization
adjusts the basal friction, and then releases them gradually. A run restarted from
a spin-up sets both to `False` to keep the optimized fields fixed.

### `[relax]`

The relaxation of the topography (`step_relax`). While active, `ytopo.topo_rel`
is `topo_rel` and its timescale ramps from `tau1` to `tau2`; after `time2`, the
`[ytopo]` values of `topo_rel` and `topo_rel_tau` apply again.

| Key | Values | |
|---|---|---|
| `topo_rel` | `ytopo.topo_rel` mode | where the ice relaxes while active (`3`: all points, `4`: the grounding line and grounding zone) |
| `tau1`, `tau2` | [yr] | relaxation timescale until `time1`, and at `time2` |
| `time1`, `time2` | [yr] | end of the `tau1` period, and end of the relaxation |
| `m` | [-] | exponent of the ramp from `tau1` to `tau2` |

### `[opt]`

The optimization (`step_optimize`), towards the observed ice thickness, of the
basal-friction coefficient `cb_ref` (`opt_cf`, between `cf_time_init` and
`cf_time_end`) and of the thermal-forcing correction `tf_corr` of the marine shelf
(`opt_tf`, between `tf_time_init` and `tf_time_end`):

| Key | Values | |
|---|---|---|
| `opt_cf` | `none`, `L21` | `cb_ref` following Lipscomb et al. (2021) |
| `opt_tf` | `none`, `L21`, `L21-points` | one `tf_corr` per basin (`tf_basins`), or one per point (`tf_sigma`, `basin_fill`) |
| `cf_init` | value | initial `cb_ref` on a cold start (`<= 0`: the till friction of the bed, `cb_tgt`) |

The method and its other parameters are described in the Yelmo docs,
[Basal friction optimization](https://fesmc.github.io/yelmo/optimization.html).

### `[output]`

Each module writes its own files, on its own grid, at the `[tm_2D]` (2D),
`[tm_2Dsm]` (small 2D) and `[tm_1D]` (time series) intervals:

| Key | Files |
|---|---|
| `write_yelmo` | `yelmo.nc`, `yelmo_sm.nc`, `yelmo_ts.nc` (and `yelmo_ts_<region>.nc` for the named regions of `[domain]`) |
| `write_isos` | `isos.nc`, `isos_ts.nc` |
| `write_shelf` | `mshlf.nc` |
| `write_surface` | the surface model's file, `<surface_method>.nc`: `chion.nc`, `smbpal.nc`, `smb_simple.nc`; none with `climate` (the climate's file has the smb, `smb_ann`, `tsrf_ann`) |
| `write_clim` | the backend's file: `snap.nc` (snapshots), `esm.nc` + `esm_ts.nc`, `rembo.nc` + `rembo_ts.nc` |
| `write_htopo` | `htopo.nc` (the hub) |
| `write_cmip`, `dt_cmip` | `yelmo_cmip.nc`, `yelmo_ts_cmip.nc`, every `dt_cmip` years (marine-shelf fields need `grid_shelf = grid_ice`) |

## Also built from this driver

`make yelmox_glaciers` previously produced a second binary from this same
`yelmox.f90` for mountain-glacier runs; it has been removed. Use `make yelmox`
with the `yelmox_Pyrenees.nml` / `yelmox_SRG.nml` configs instead.
