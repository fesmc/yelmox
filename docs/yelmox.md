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
- **Configs:** one base par file per domain and climate backend,
  `yelmox/yelmox_<base>.nml`: Greenland, GreenlandPaleo, Antarctica, North, LIS,
  Pyrenees, SRG, esm_Greenland and esm_Antarctica ([ESM forcing](climate-esm.md)),
  rembo_Greenland ([REMBO](climate-rembo.md)). An experiment is a `-p` list on a
  base, in its scripts in `scripts/<exp>/` (arrays shared by several scripts in
  `params.sh`); see [Forcing library](#forcing-library).
- **Experiments:** [Greenland](experiment-greenland.md),
  [Antarctica](experiment-antarctica.md), [ISMIP7 projections](experiment-ismip7.md),
  [1pctCO2](experiment-1pctco2.md), [TIPMIP](experiment-tipmip.md),
  [Antarctic paleo](experiment-ant-paleo.md), [Greenland paleo](experiment-grl-paleo.md),
  [Laurentide](experiment-laurentide.md), [Mountain glaciers](experiment-mountain.md).

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
| `snapshots` | atmosphere + ocean profiles: a reference snapshot plus the terms of `[forcing]` (coefficient × snapshots or uniform anomalies) | [Snapshots](climate-snapshots.md) |
| `esm` | reference climatology + Earth-system-model anomalies; atmosphere or smb, ocean at the shelf base, subglacial discharge | [ESM forcing](climate-esm.md) |
| `rembo` | REMBOv1 atmosphere + smb; ocean from snapshots (ocean-only snapshot) | [REMBO climate](climate-rembo.md) |

## Stepping order

The driver owns the timeline (`ts`) and the shared sea level (`bsl`), and advances
the domain once per step with the coupling sequence written out in the time loop:

```fortran
call step_relax(dom, ts)          ! topography relaxation (relax)
call step_optimize(dom, ts, dtt)  ! cb_ref/tf_corr optimization (opt)
call step_isostasy(dom, ts, bsl)  ! bedrock + sea level (dt_isos cadence)
call couple_to_yelmo(dom)         ! bedrock now; smb + shelf melt lag one step
call step_ice_bnd(dom, ts, coefs) ! prescribed ice in the fixed/relaxation cells ([ice_bnd])
call step_icesheet(dom, ts)       ! yelmo_update
call couple_yelmo_to_htopo(dom)   ! hi-res geometry from the models
call step_climate(dom, ts, dtt, coefs) ! climate (dt_clim cadence)
call step_surface(dom, ts)        ! surface mass balance + temperature
call step_shelf(dom, ts)          ! shelf-base melt + temperature
```

`couple_to_yelmo` assembles the Yelmo boundary state: isostasy from this step,
smb and shelf melt from the **previous** step (a one-step coupling lag);
`step_climate`, `step_surface` and `step_shelf` then produce the forcing
consumed on the next step.
The climate is refreshed on the `comps.dt_clim` cadence, isostasy on the
`comps.dt_isos` cadence (the bedrock is held in between); the smb every step.

## Forcing: coefficients and terms (`[forcing]`, `[coef_<name>]`, `[term_<name>]`) {#forcing}

The forcing of a run is one equation set, and the par file is its
transcription (design: `docs/_forcing-terms.md`). At the **forcing time**
τ = `time_rel` (set by `[ctrl] tstep_method`: `const` holds `tstep_const`,
`rel` is the model time, `cal` the calendar), the climate is a **reference**
snapshot plus **terms**, each a **coefficient** (a scalar of τ) times a **pattern**:

```
T^sl(x,m,τ) = T_0^sl + Σ_k α_k(τ) [T_k1^sl − T_k0^sl]          + Σ_u f_ta,u α_u(τ)
P^sl(x,m,τ) = P_0^sl · Π_k [1 + α_k(τ) (P_k1^sl/P_k0^sl − 1)]  · Π_u exp(β_p f_ta,u α_u(τ))
Θ(x,z,τ)    = Θ_0 (or the coefficient of a replace term) + Σ_k α_k [Θ_k1 − Θ_k0] + Σ_u f_to,u α_u
S(x,z,τ)    = S_0 (or ...)                          + Σ_k α_k [S_k1 − S_k0] + Σ_u f_so,u α_u
```

at sea level (the snapshots reduced with the lapse rates Γ and β_p of
`[snapshots]`), then moved to the model surface; k are the terms with a
`snapshots` pattern, u those with a `uniform` pattern. The rules are fixed by
field: `tas`, `to`, `so` add, `pr` multiplies (the ratio of the snapshots; in a
uniform term exp(β_p ΔT)). A term without `snap0` adds `snap1` itself (a given
anomaly or ratio). The same equations hold for every backend that takes terms
(snapshots; rembo takes the uniform `tas` terms as its summer anomaly; esm
takes ocean terms at the shelf base). A **parameter** is a fixed number of the
par file: the gains `f_ta`, `f_to`, `f_so` of a uniform term, β_p and the lapse
rates of `[snapshots]`.

```fortran
&forcing
    reference = "pd"                 ! the reference snapshot [snapshot_pd]
    coefs     = "alpha" "dTo"        ! [coef_alpha], [coef_dTo]
    terms     = "lgm" "ocean"        ! [term_lgm], [term_ocean]
    var_defs  = "input/clim_db.nml" "input/forcing_db.nml"   ! the default
/

&coef_alpha                          ! a tsgen coefficient (unset keys: template [coef])
    method      = "series"
    series_file = "input/alpha_combined_125kyr_interp.dat"
/

&coef_dTo                            ! derived: set by the climate each update
    method = "derived"
    from   = "air_anomaly"           ! factor x mean sea-level air anomaly over region
    factor = 0.25
    region = "all"
    units  = "K"
/

&term_lgm                            ! tas + alpha (lgm - picontrol), pr (1 + alpha (lgm/picontrol - 1))
    coef   = "alpha"
    snap0  = "picontrol"
    snap1  = "lgm"
    fields = "tas" "pr"
/

&term_ocean                          ! to + dTo at every depth
    coef    = "dTo"
    pattern = "uniform"
    fields  = "to"
/
```

This reads T = T_pd + α (T_lgm − T_pic), P = P_pd [1 + α (P_lgm/P_pic − 1)],
Θ = Θ_pd + 0.25 ⟨T^sl − T^sl_pd⟩ (Greenland paleo A). The coefficients are named
by their symbol in the equation (`alpha`, `beta`, `gamma`, ...), the record
each loads is the setup's: `[coef_alpha]` is `alpha_combined_125kyr` in
`yelmox_GreenlandPaleo.nml` and `alpha_orbital_r2020` (normalised) in
`yelmox_Greenland.nml`, so a term on `alpha` (`term_lgm`) is the same in both.
The library `input/forcing_db.nml` holds shared groups ([Forcing
library](#forcing-library)), so a par file lists names; an experiment adds or
removes terms with `-p forcing.terms==lgm,ocean`.

**Coefficients** are named scalars on the forcing time, **driver-owned**: the
program holds one registry (`coefs_class`, `coefs`, libs/kryos_forcing.f90),
advances every coefficient once per step and passes it to `step_climate` (and the
cold start); a multi-domain driver shares one registry. Each is a `tsgen`
time-series generator (`method`: `const`, `ramp-slope`, `ramp-time`,
`ramp-time-step`, `sin`, `series` from a table or netCDF file, optionally
normalised with `norm_t0`/`norm_t1` to 0 at t0 and 1 at t1, each one time or a
window t0 t1 whose mean is taken, `exp`, the PI/PID feedback controllers on the
response, the ice volume in Gt; see the fesm-utils tsgen docs), **scaled**
(`method = "scaled"`: `factor` × the tsgen coefficient `coef`, set right
after it each step, so a term on the atmosphere may use it; e.g. a regional
amplification of a global anomaly, `[coef_dTa_greenland]` of the library), or
**derived** (`method = "derived"`), computed by the model after the atmosphere
of the update is composed:

| `from` | value | keys |
|---|---|---|
| `air_anomaly` | `factor` × the mean over `region` of the sea-level annual air temperature anomaly of this update (no lag; without the lapse-rate cooling of the ice sheet) | `factor`, `region` |
| `obm_north`, `obm_south` | the box temperature of the bipolar ocean box model, minus its initial value with `anomaly = True` ([yelmox_bipolar](yelmox-bipolar.md)) | `anomaly` |

A term on a derived coefficient may touch `to` and `so` only (an error at init
otherwise); `air_anomaly` needs a backend with a sea-level air temperature
(snapshots). `region` is `"all"` or a selection expression of the domain's
regions on the climate grid (`[regions]`, e.g. `"zone:open_ocean"`,
`"~zone:land"`; see the fesm-utils regions docs); its cell count is printed at
init. The run stops when a coefficient trips its kill switch (`with_kill`), and
writes a restart bundle `restart-<n>` each time a coefficient moves by more than
its `restart_every_df`. Each tsgen coefficient restarts from
`coef_<name>_restart.nc` in the bundle (a scaled or derived one is recomputed; an older
bundle's `series_<name>_restart.nc` is still read, with a note); every
coefficient is written to the run's 1D file `forcing_ts.nc` (`[tm_1D]`) under
its name (tsgen ones also `<name>_df_dt`), labelled with `units`, its `kind`
and source as attributes, next to the barystatic sea level `bsl`. The names
`time` and `bsl` are therefore not coefficient names.

**Terms** are one coefficient times one pattern on an explicit list of `fields`:

| `pattern` | adds | keys |
|---|---|---|
| `snapshots` (default) | α (snap1 − snap0), or α snap1 with `snap0 = ""` (a sea-level anomaly or ratio, e.g. a reconstruction in time) | `snap0`, `snap1`, `fields` |
| `uniform` | `f_ta` α to `tas` (and exp(β_p `f_ta` α) to `pr`), `f_to` α to `to`, `f_so` α to `so` | `fields`, `f_ta`, `f_to`, `f_so` (default 1) |
| `replace` | the field = α at every depth, in place of the reference; the other terms add to it | `fields` (`to`, `so`), one per field |

`coef = ""` is the constant 1 (the full anomaly of a reconstruction). The
snapshots of a run are the reference and the `snap0`/`snap1` of its terms: no
list. Checked at init, with an error naming the term: every coefficient and
snapshot named exists, every field is bound by its snapshots and the reference,
a `replace` once per field, a derived coefficient on the ocean only, every
coefficient of `[forcing] coefs` used (by a term, `[ice_bnd]`, the OBM forcing
of yelmox_bipolar, or as the `coef` of a scaled coefficient that is); a backend that
supplies the smb (esm with `surface_method = "climate"`) refuses terms on
`tas`/`pr`. The value of each term's coefficient is written to the climate's 2D
file as `term_<name>`.

**Order** of a climate update: the backend composes the atmosphere with the
snapshot terms, then the uniform terms on `tas`/`pr`; the derived coefficients are
set (`air_anomaly` from the composed sea-level temperature); then the ocean
terms (`replace`, then snapshot and uniform); then the transforms `[taper]`
and `[bias]`; then the move to the surface grid's elevation. The ocean terms
are re-applied by `yelmox_bipolar` after the box model sets its coefficients.

Every named group (`coef_`, `term_`, `snapshot_` and the varslice groups of
the fields, the taper weight, the bias fields, the ice history) is read from
the par file if it has it, else from the first `[forcing] var_defs` file that
does; unset keys come from the templates of `input/yelmox_defaults.nml`
([Configuration](#configuration)).

## Regional taper to the reference (`[taper]`)

The climate transform after the terms and before `[bias]`: the atmosphere
towards the reference climate with a static weight field w in [0, 1] on the
climate grid,

    now = ref + w*(now - ref)

for the monthly air temperature and precipitation, their annual means and the
sea-level air temperature (those the backend supplies); w = 1 keeps the
climate, w = 0 is the reference (e.g. outside the coverage of a
reconstruction). The reference is taken at the current
surface: a reference at its own elevation (snapshots: the surface of the
reference snapshot) is moved there with the backend's monthly lapse rate and
`beta_p`, as the climate itself. The ocean is unchanged. Refused for a backend
that supplies the smb itself (esm with its smb; rembo).

`outside` sets the climate where w < 1:

- `"ref"` (default): towards the reference, as above;
- `"mean"`: towards the reference plus the mean anomaly of the cells with
  w = 1, per month,

      tas = tas_ref + w*(tas - tas_ref) + (1 - w)*dT(m)
      pr  = w*pr + (1 - w)*pr_ref*f(m)

  with dT(m) the mean temperature difference and f(m) the ratio of the
  precipitation totals over the cells with w = 1 (both against the reference
  at the current surface); `ta_ann` and `tsl_ann` take the mean of dT(m),
  `pr_ann` the ratio of the annual totals. A weight without cells at w = 1
  stops the run.

```
&taper
    active   = True
    weight   = "taper_w"  ! group of w [1]: one 2D field (time = 0 0 0 0 1); the par file, then [forcing] var_defs
    outside  = "ref"      ! ref | mean: where w < 1, towards the reference (plus the mean anomaly of the cells with w = 1)
/
```

The weight is a varslice group on the climate grid or remapped online
(`remap`, `grid_src`); values outside [0, 1] (or missing) stop the run. It is
written to the climate's 2D file (`taper_w`). `input/clim_db.nml` holds
`taper_b18_b20`, the weight of the Greenland reconstructions on GreenlandPaleo
(1 in their coverage, 0 beyond 100 km; `outside = "mean"` in `grlpal_recon`
of `scripts/grl-paleo/params.sh`).

## Climate bias correction (`[bias]`)

A spatial correction of the atmosphere, the last climate transform, applied
on the climate grid to the monthly air temperature and precipitation, before
the move to the surface grid's elevation; refused for a backend that supplies
the smb itself (esm with its smb; rembo):

    tas = tas + f_tas*dtas,    pr = pr*(1 + f_pr*(fpr - 1))

The annual means follow (temperature additively, precipitation by the ratio
of the annual totals); the reference climate and the ocean are unchanged. The
group is optional; without it, or with `active = False`, nothing changes.

```
&bias
    active   = True
    dtas     = "bias_dtas"    ! group: temperature correction [K]; "" = none (the par file, then [forcing] var_defs)
    fpr      = "bias_fpr"     ! group: precipitation factor [1]; "" = none
    f_tas    = 1.0            ! weight of dtas
    f_pr     = 1.0            ! weight of fpr
/
```

`dtas` and `fpr` are varslice groups: one field for all months
(`time = 0 0 0 0 1`) or 12 monthly fields (`time = 0 0 0 0 12`), on the
climate grid or remapped online (`remap`, `grid_src`). A group with a time axis
(`time = 1 x0 x1 dx nsub`) is a time-dependent bias: it is sliced at the
forcing time at every climate update, as a time-varying snapshot. The weighted
annual fields are written to the climate's 2D file (`bias_dtas_ann`,
`bias_fpr_ann`; with a time axis also `bias_dtas_ann_now`, the slice applied).

## Insolation and orbit (`[orbit]`)

The surface models that take the top-of-atmosphere insolation (chion's ITM and
BESSI, smbpal's ITM) compute it with the [insol](https://github.com/fesmc/insol)
library (a dependency, `libinsol.a`), the daily mean in closed form from the
orbital elements of Laskar et al. (2004). `[orbit]` (one per domain,
`[orbit_north]` in yelmox_bipolar) gives the elements; both surface models read
it:

```
&orbit
    method = "fixed_time"   ! forcing_time: at the forcing time time_rel; fixed_time: at time; fixed_orbit: ecc, obl, per
    time   = 0.0            ! [yr] model time of the orbit (fixed_time; relative to 1950, past < 0)
    ecc    = 0.016702       ! [1]   eccentricity (fixed_orbit)
    obl    = 23.4393        ! [deg] obliquity (fixed_orbit)
    per    = 102.918        ! [deg] longitude of perihelion (fixed_orbit)
    fldr   = "input"  ! folder of the orbital tables
/
```

`forcing_time` follows the forcing time (`ts%time_rel`: transient insolation
with `tstep_method = "rel"`); `fixed_time` holds the orbit of one time
(present day by default, the 1950 orbit); `fixed_orbit` takes the given
elements (the defaults are the tables' at 1950). A fixed orbit is the one
exception to the single forcing time. Each model recomputes its insolation
only when the elements change, so once per run unless
`method = "forcing_time"`. The par files list `method` and `time`; the other
keys come from `input/yelmox_defaults.nml`.

## Prescribed ice (`[ice_bnd]`)

Yelmo prescribes the ice thickness in its fixed cells (`mask_ice` = 1:
`H_ice = H_ice_ref`) and relaxes it towards `H_ice_ref` in its relaxation cells
(`tau_relax > 0`, with `ytopo.topo_rel = -1`), set by `[yelmo_masks]`: a strip
at the edge of the region where ice is allowed (`border_width = N M`: N fixed
cells, then M relaxation cells with `relax_tau`), the expressions
`mask_ice_fixed` and `relax`, or `yelmo.mask_border`. `[ice_bnd]` writes the
thickness of an ice history into `H_ice_ref` in these cells, each step before
`yelmo_update` (`step_ice_bnd`) and at the cold start; the interior keeps its
`H_ice_ref` (the present-day ice):

    file:    H_ice_ref = H_ice(time_rel)        (interpolated; the first or last slice outside the file)
    scaled:  H_ice_ref = f * H_ice(slice_time)  (f a coefficient, e.g. before the file starts)

```
&ice_bnd
    active     = True
    method     = "file"               ! file | scaled
    H_ice      = "paleomist_H_ice"    ! group of the ice thickness [m] (the par file, then [forcing] var_defs)
    coef       = ""                   ! scaled: the coefficient f ([forcing] coefs)
    slice_time = -21e3                ! [yr] scaled: time of the slice
/
```

`input/clim_db.nml` holds `paleomist_H_ice` (PaleoMIST 1.0, 80-0 ka every
2.5 kyr) and `glac1d_H_ice` (GLAC-1D, 26-0 ka every 100 yr), the IceHistory
products on the ice grid. The forcing time is `time_rel` (a constant-time
run holds the ice of `tstep_const`). At the cold start the isostasy
reference is taken before (the present-day ice: the prescribed ice of another
time is a load). With a time-varying border use `[barysealevel] method =
"file"`: under `fastiso` the prescribed ice would add to the barystatic sea
level. `H_ice_ref` is written to `yelmo.nc`. Example: `yelmox_GreenlandPaleo.nml`
(GRL-PAL, `border_width = 2 4`, `mask_ice_dynamic = "all"`).

## Configuration

A run is set by these namelist groups. Every key is required (a missing key stops
the run), except the keys of an unselected `init_method` and the `[relax]` and
`[opt]` groups when their switch in `[sim]` is off.

The forcing groups have defaults: `input/yelmox_defaults.nml` declares every
key of `[forcing]`, `[coef_<name>]` (template `[coef]`: the tsgen keys,
the derived keys, `units`, `restart_every_df`), `[term_<name>]` (template
`[term]`), the snapshot backend (`[snapshots]`, template `[snapshot]` for
`[snapshot_<name>]`), the climate transforms (`[taper]`, `[bias]`), the orbit
of the insolation (`[orbit]`), `[ice_bnd]` and `[obm_forcing]`
(yelmox_bipolar), with the domain suffix of a multi-domain driver
(`[forcing_north]` to `[forcing]`). A par file lists only what differs; an
absent group takes every default (a transform is then inactive). Each read goes through
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
| `[domain]` | the domain: name, physical constants, grid of every component, basin set of the hub (`basins`) (see [Multigrid coupling](multigrid.md)) |
| `[topo]`, `[regions]` | the hub topography; the regions, zones, basins and named masks (`masks`, `mask_<name>`) |
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
    write_lossless  = False             ! all output lossless (ignore per-field precision)
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

### Forcing library

`input/forcing_db.nml` holds the named forcing groups that several par files
share: snapshots `[snapshot_<name>]` (PMIP3 `picontrol`, `lgm`; the Greenland
reconstructions `b18_base`, `b18`, `b20`), coefficients `[coef_<name>]`
(`alpha`, the default record, read when the par file has none; `dTo`,
`ramp`, `hyst`, `dT`, `obm_north`, `obm_south`), terms `[term_<name>]`
(`lgm`, `b18`, `b20`, `ocean`, `ramp`, `hyst`, `dT`, `obm_north`, `obm_south`) and
elsa layer settings (`elsa_grl_lgp`, `elsa_grl_deglac`, `elsa_ant_paleo`).
`[forcing]` (the per-run selection), `[snapshots]` (the backend's rules) and
`[snapshot_pd]` (domain-specific) stay in the par file.

A named group is read from the par file if it has it, else from the first
`[forcing] var_defs` file that does (default `"input/clim_db.nml"
"input/forcing_db.nml"`: the fields, then the library). elsa reads its group
from the file `[ytrc] elsa_nml` names (`input/forcing_db.nml`), selected by
`elsa_group`. A listed group found nowhere stops the run.

To change a library group for one run, define it in the par file (it wins; LIS
and North redefine `coef_alpha` as the constant 1), or set one of its keys
with `-p` (e.g. `-p coef_alpha.f_min=0`): runme copies the whole group into
the staged par file, then applies the override (`par_library` in
`.runme/info.json`). `par_defaults` (`input/yelmo_defaults.nml`,
`input/yelmox_defaults.nml`) lets `-p` set a key of a group the par file lacks,
e.g. `-p bias.active=True` or `-p ytrc.use_elsa=True`. `-p` on a library group
needs runme 0.5.26 or later.

Scripts follow one rule: the settings that define the experiment go in
commented `-p` arrays (in `scripts/<exp>/params.sh` when several scripts share
them); run control goes in the invocation (run step, times, restart, grid,
output intervals, queue).

The [forcing builder](forcing-builder.md) builds the forcing in the order of
its equation (time, reference and terms, the records and fields that populate
them, the transforms), checks that every piece is there, and writes the run
out as the par-file groups and as the `-p` list with the `runme` command
against a base; `scripts/forcing-builder/explain.jl` resolves a base plus
`-p` the same way on the command line.

### `[ctrl]`

| Key | Values | |
|---|---|---|
| `run_step` | group name | the group holding the timeline: `"ctrl"` for `[ctrl]` itself, or a phase group such as `[spinup]` or `[transient]`, so one par file can hold several run phases |
| `calendar`, `calendar_ref` | bool, [yr CE] | timeline in calendar years, against the reference year `calendar_ref` (ESM runs) |
| `restart` | folder, `"None"` | the restart bundle to start from; `"None"` = cold start |
| `write_lossless` | bool | optional (default `False`): write all output lossless, ignoring the per-field precision (e.g. for bit-for-bit comparisons of runs) |

The timeline group holds `tstep_method` (`"const"`: the time advances and the
forcing time stays at `tstep_const`; `"rel"`: the forcing follows the time, e.g.
for paleo records; `"cal"`: calendar time), `tstep_const`, `time_init`,
`time_end` and the main time step `dtt` [yr].

A restart bundle is a folder `restart-<kyr>-kyr/`, written on the `[tm_rst]`
schedule and always at `time_end`. It holds one restart file per stateful model,
the shared sea level and the transient-forcing state (in `yelmox_bipolar` also the
ocean box model, with each domain in a subfolder named after it). Pass `restart`
as an absolute path: the run starts inside its own folder. A restart continues the
run bit for bit: the state is restored, not stepped again (the surface restart
files carry the annual fields of the last cycle); only the climate is re-derived
from it, which equals the run's own when the restart time lies on its `dt_clim`
cadence (counted from the start of the run, as are `dt_isos` and the `[relax]`,
`[opt]` time windows). With `dt_clim <= 0` the climate is re-derived once at the
restart, from the restored geometry, so it is not the cold-start climate of the
run that wrote the bundle.

### `[comps]`

The components: which are active, with which model, how often.

| Key | Values | |
|---|---|---|
| `with_ice_sheet`, `with_isostasy`, `with_climate`, `with_surface`, `with_shelf` | bool | components in the coupling sequence: the ice sheet, isostasy, the climate (atmosphere + ocean), the surface (mass balance + temperature) and the shelf base (melt + temperature); `with_surface` and `with_shelf` need `with_climate` |
| `climate` | `snapshots`, `esm`, `rembo` | the climate backend |
| `surface_method` | `chion` (default), `smbpal`, `smb_simple`, `climate` | chion or smbpal (ITM or PDD), both from the climate's temperature and precipitation; smbpal's ITM spins its snowpack up for `[smbpal] time_equil` years at the cold start; chion runs the ITM, BESSI or PDD snowpack (`[chion] model`) in daily steps (`[surface_chion]`, `[chion]`; ITM's parameters are the group `[chion] nml_itm` names: smbpal's `[itm]`, or for Greenland `[itm_chion]` (`alb_ice = 0.31`, calibrated to MAR); BESSI gets the surface shortwave as TOA insolation times `trans_sw`, the TOA insolation itself (chion's cloud-proxy longwave), the ice thickness, a constant `wind_speed` and `rel_hum`, and the air pressure from the surface elevation; its constants are `[chion_const]`, for Greenland `alpha_ice = 0.55`, calibrated to MAR; BESSI runs pass `surface_chion.time_equil=20`, converged within 0.3 Gt/yr of the par files' 100 yr, which ITM needs; PDD's parameters are the group `[chion] nml_pdd` names, `[pdd_chion]` with `pdd_method = "cg05"` (Calov and Greve, 2005), smbpal's PDD; the surface temperature of ITM and PDD is chion's `chion_tsrf_annual`, the annual-mean air temperature warmed by `firn_fac` times the annual refreezing, capped at freezing over ice); smb_simple (needs a sea-level air temperature: snapshots); the climate's own smb (esm, rembo; required by rembo) |
| `dt_clim` | [yr] | climate update interval; `<= 0`: updated only at the cold start (and re-derived once on a restart) |
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
| `use_negis` | bool | NEGIS basal-friction modification (`[negis]`): `cb_ref` is set each step (Yelmo's `ytill.method = -1`) as the till friction of the bed with the NEGIS parts scaled, the centre by `cf_centre`, the north part by `cf_north` before 4 ka BP, the south part by `cf_south` from 4 ka BP, all parts by `cf_0` before 11 ka BP ramping linearly to `cf_1` at 0 ka, on the forcing time; not with the `cb_ref` optimization (`opt_cf`) |
| `scale_glacial_smb` | bool | reduce negative glacial smb (`[glacial_smb]`: `dt_lgm`, `lat_lim`, `fac_lim`), by a glacial index from the climate's annual air temperature against its reference; refused with `climate = rembo` (no reference) |
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
(`opt_tf`, between `tf_time_init` and `tf_time_end`). `cb_ref` is optimized on
`grid_ice`, `tf_corr` on `grid_shelf`, where it lives: the Yelmo-grid inputs
(current and observed ice thickness and grounding, `dHidt`) are remapped onto
`grid_shelf`, the basins (and the ids of `tf_basins`) are those of
`[marine_shelf] basins`, and `tf_corr` is never remapped.

| Key | Values | |
|---|---|---|
| `opt_cf` | `none`, `L21` | `cb_ref` following Lipscomb et al. (2021) |
| `opt_tf` | `none`, `L21`, `L21-points` | one `tf_corr` per basin (`tf_basins`), or one per point (`tf_sigma`, `basin_fill`) |
| `cf_init` | value | initial `cb_ref` on a cold start (`<= 0`: the till friction of the bed, `cb_tgt`) |

The method and its other parameters are described in the Yelmo docs,
[Basal friction optimization](https://fesmc.github.io/yelmo/optimization.html).

### `[output]`

Each module writes its own files, on its own grid, at the `[tm_2D]` (2D),
`[tm_2Dsm]` (small 2D) and `[tm_1D]` (time series) intervals. Every key is
required. The driver's forcing, shared by its domains (the coefficients and
the barystatic sea level), goes to `forcing_ts.nc` in the run folder on the
`[tm_1D]` cadence, whatever the switches. A component switched off in `[comps]` writes no files whatever its
switch (`write_isos`, `write_shelf`, `write_surface`, `write_clim` need
`with_isostasy`, `with_shelf`, `with_surface`, `with_climate`):

| Key | Files |
|---|---|
| `write_yelmo` | `yelmo.nc`, `yelmo_sm.nc`, `yelmo_ts.nc` (and `yelmo_ts_<name>.nc` for the named masks of `[regions]`) |
| `write_isos` | `isos.nc`, `isos_ts.nc` |
| `write_shelf` | `mshlf.nc` |
| `write_surface` | the surface model's file, `<surface_method>.nc`: `chion.nc`, `smbpal.nc`, `smb_simple.nc`; none with `climate` (the climate's file has the smb, `smb_ann`, `tsrf_ann`) |
| `write_clim` | the backend's file: `snap.nc` (snapshots), `esm.nc` + `esm_ts.nc`, `rembo.nc` + `rembo_ts.nc` |
| `write_htopo` | `htopo.nc` (the hub) |
| `write_cmip`, `dt_cmip` | `yelmo_cmip.nc`, `yelmo_ts_cmip.nc`, every `dt_cmip` years (marine-shelf fields need `grid_shelf = grid_ice`) |

The 2D fields are written at a fixed precision (ncio bit rounding, attributes
`precision_abs`/`precision_rel`), so that they compress well: 1e-3 m for
elevations, 1e-3 K for temperatures, 1e-3 m/yr (or 1e-4 relative) for
velocities, 1e-4 m/yr for mass balance rates, 1e-4 for fractions, 1e-4 relative
for stresses, precipitation and the other fields; masks, time series, restarts
and the static fields are lossless. `[ctrl] write_lossless = True` writes
everything lossless.

## Also built from this driver

`make yelmox_glaciers` previously produced a second binary from this same
`yelmox.f90` for mountain-glacier runs; it has been removed. Use `make yelmox`
with the `yelmox_Pyrenees.nml` / `yelmox_SRG.nml` configs instead.
