---
title: "Snapshots"
aliases:
  - climate-snap.html
---

The snapshot backend `snapshots` builds the climate from a reference state plus
climate **snapshots** (e.g. present day, pre-industrial, LGM), blended in time
by scalar **indices** (e.g. a glacial index), as one blend model over any number
of snapshots, fields and axes, on the fesm-utils `varslice` reader and `tsgen`
series. It supplies the atmosphere (monthly `tas`, `pr`, lapse-rate corrected to
the current surface) and the ocean as depth profiles (`to_ann`, `so_ann`), which
the marine shelf interpolates to the shelf base. It runs on `grid_clim`, reads
the `[snap]` groups and also supplies the ocean of the rembo climate (an
ocean-only reference snapshot).

- **Program:** `yelmox` (`make yelmox`), backend in `libs/snapshots.f90`, behind
  `libs/yelmox_climate.f90`.
- **Select:** `[comps] climate = "snapshots"`.
- **Output:** with `[output] write_clim`, `snap.nc`, the climate as applied
  (after the climate modifiers): the annual air temperature and its summer
  mean (`t2m_ann`, `t2m_sum`), the annual precipitation, the ocean profiles
  (`to_ann`, `so_ann` on `depth`), the position on each axis (`axis_<name>`),
  the surface elevation, the monthly atmosphere (`tas`, `pr`), `tsl_ann` and
  the diagnostics of the active climate modifiers.

The climate modifiers `[anom]`, `[ocn_fraction]`, `[taper]` and `[bias]` (see
[yelmox](yelmox.md)) apply after the backend, in this order; `[ocn_obm]` is
applied by `yelmox_bipolar`.

## Configuration

Example: `yelmox/yelmox_Greenland.nml`, present day and PMIP3 piControl and LGM
on one axis; its ocean anomaly comes from `[ocn_fraction]`.

```fortran
&snap
    var_defs  = "input/clim_db.nml"   ! databases of the input files (default)
    ref_name  = "pd"                  ! the reference snapshot (default)
    snapshots = "pd" "picontrol" "lgm"
    axes      = "at"
    lapse     = 0.0080 0.0065         ! [K m-1] annual, summer (default)
    f_p       = 0.06                  ! [K-1] precipitation scaling (default 0.05)
/
```

Every key of the `[snap*]` groups has a default in `input/yelmox_defaults.nml`
(see [Defaults](#defaults)): `combine = "anomaly"` (snapshots as anomalies from
the reference), `fields = "tas" "pr" "to" "so" "zs"`, `grid_global =
"GLOBAL-1DEG"`, ... Three kinds of group follow:

- **Fields**, `[snap_field_<name>]`: `blend` (`linear`, `ratio`, `const`;
  default `linear`) and the `axes`, one or more names (default `""` = none).
  The field name (`tas`, `pr`, `to`, `so`, `zs`) decides what the field is. A
  field without axes follows the reference: its group can be left out.
- **Snapshots**, `[snap_<name>]`: `exp` (the experiment, `{exp}` in the
  database paths), `monthly` (default `True`), `clim_sigma`/`ocn_sigma` ([km]
  Gaussian smoothing of the atmosphere and the ocean at loading, the land mask
  following the smoothed elevation; default 0), `tas_season` and, per field,
  the database group(s) holding it (see [Database](#database)).
- **Axes**, `[snap_axis_<name>]`: two snapshots `lo` -> `hi` and the `index`
  giving the position on the axis: 0 at `lo`, 1 at `hi`. The index is a
  driver-owned `[idx_<name>]` (a `tsgen` generator of any method: `const`,
  ramps, `series`, ...; see [yelmox](yelmox.md)), listed in `[indices]`.

### Database

The snapshots bind varslice groups (see [Database namelists](database-namelists.md))
by name, from the par file itself if it has the group (searched implicitly
first, so it never has to name itself and its groups win), else from the
files listed in `[snap] var_defs`, searched in order (first file that has the
group). The snapshot groups `[snap_<name>]` and the axis groups
`[snap_axis_<name>]` are looked up the same way, and the index groups
`[idx_<name>]` through `[indices] var_defs` (default none). The shared
database `input/clim_db.nml` holds one group per product field; a par file can
add run-specific groups:

```fortran
&snap
    var_defs = "input/clim_db.nml"   ! the default
/
&my_tas                              ! in the par file: found before clim_db.nml
    filename = "ice_data/..."
    ...
/
```

A group missing from the par file and every var_defs file stops the run,
naming the group and the files.
The paths of the database take, besides `{domain}` and `{grid_name}`, two
keywords:

- `{grid_global}`: the global grid of the paleoclimate products, `[snap]
  grid_global` (default `GLOBAL-1DEG`);
- `{exp}`: the experiment of the snapshot, its `exp`, so one group serves every
  experiment of a product family.

Present-day references are regional (`{domain}/{grid_name}`, the resolution of
the climate grid): MAR v3.14 1961-1990 (`mar314_1961_1990_{tas,sf,rf,pr,zs}`),
RACMO2.3p2 1981-2010 (`racmo23p2_1981_2010_{tas,tsrf,pr,sf,zs}`), ERA5
1981-2010 (`era5_1981_2010_{tas,pr,sf,zs}`), ORAS5 1981-2010
(`oras5_1981_2010_{to,so}`, monthly climatology) and the ISMIP7 ocean
(`ismip7_{to,so}`). The paleoclimate snapshots are the global products,
remapped online onto the climate grid (`remap = "bilin"`, the map cached in
`maps/`):

| Product (FesmData) | `{exp}` | Groups |
|---|---|---|
| `Paleoclimate-PMIP3-{exp}-mean` | lgm, midHolocene, piControl | `pmip3_mean_{tas_ann,tas_jja,tas_djf,pr_ann,zs}` |
| `Paleoclimate-PMIP4-{exp}-mean` | lgm, midHolocene, piControl | `pmip4_mean_{tas_mon,pr_mon,tas_ann,tas_jja,tas_djf,pr_ann,zs}` |
| `Paleoclimate-CLIMBER3a-{exp}`, `-ocean` | present, lgm_1p7weak, lgm_1p7strong | `c3a_{tas_ann,tas_jja,pr_ann,zs}`, `c3a_to` |
| `Paleoclimate-CLIMBERX-Willeit2025-{exp}`, `-ocean` | piControl, interstadial, stadial, heinrich | `cx_{tas_mon,pr_mon,tas_ann,pr_ann,zs}`, `cx_{tas_mon_anom,pr_mon_anom}` (anomalies to piControl, for `lo = ""` axes; `pr_mon_anom` is a difference in mm d-1, not a ratio), `cx_{to,so}` |

The Greenland reconstructions are regional products on GRL-16KM, remapped
online (`grid_src = "GRL-16KM"`) and filled from the nearest cell outside
their box (`fill = "nearest"`, e.g. on GRL-PAL):

| Product (FesmData) | Time axis | Groups |
|---|---|---|
| `Paleoclimate-Buizert2018` | monthly, every 10 yr, -22005 to 55 yr (22005 BP to 2005 CE) | `b18_{tas,pr}` (time-varying), `b18_zs` (static) |
| `Paleoclimate-Buizert2018-1981-2010` | none (its 1981-2010 baseline) | `b18_1981_2010_{tas,pr,zs}` |
| `Paleoclimate-Badgeley2020` | annual, every 50 yr, 20 ka to 0 | `b20_tas_anom` [K], `b20_pr_frac`, `b20_pr_frac_low`, `b20_pr_frac_high` (anomaly and ratios relative to 1850-2000 CE, for `lo = ""` axes) |

The database also holds the groups of the climate modifiers (`[bias]`:
`bias_mar_era5_1981_2010_{dtas,fpr}`; `[taper]`: `taper_b18_b20`) and the ice
histories of `[ice_bnd]` (`paleomist_H_ice`, `glac1d_H_ice`).

**Time-varying groups.** A group with a time axis (`time = 1 x0 x1 dx nsub`)
is sliced at the forcing clock `time_rel` at the start and at every climate
update: linear interpolation between the two bracketing records, the first or
last record outside the axis (varslice `extrap`). Its snapshot is then reduced
to sea level again (and the reference reset, if it is the reference). A group
without a time axis (`time = 0 ...`) is read once at the start.

```fortran
&snap_pd
    tas = "mar314_1961_1990_tas"
    pr  = "mar314_1961_1990_sf" "mar314_1961_1990_rf"
    zs  = "mar314_1961_1990_zs"
    to  = "oras5_1981_2010_to"
    so  = "oras5_1981_2010_so"
/

&snap_lgm
    exp        = "lgm"
    monthly    = False
    clim_sigma = 50.0
    tas = "pmip3_mean_tas_ann" "pmip3_mean_tas_jja"
    pr  = "pmip3_mean_pr_ann"
    zs  = "pmip3_mean_zs"
/
```

### Defaults

`input/yelmox_defaults.nml` declares every key of the `[snap*]` groups, the
indices and the climate modifiers, with its default. The named groups fall
back to a template group: `[snap_field_<name>]` to `[snap_field]`,
`[snap_<name>]` (a snapshot) to `[snap_snapshot]`, `[snap_axis_<name>]` to
`[snap_axis]`, `[snap]` to `[snap]`, with the domain suffix of a multi-domain
driver (`[snap_north_field_tas]` to `[snap_field]`). A par file lists only what
differs; a key the defaults file does not declare stops the run
(`nml_validate`), and `scripts/check_par_nml.py` reports it offline.

### Fields and blending

**Fields are optional per snapshot.** A snapshot binds the fields it names
(a key absent or `""` binds nothing). The backend builds what the reference
snapshot (`ref_name`) binds: the atmosphere (`tas` and `pr`, with the
elevation `zs`) and/or the ocean (`to`, optionally `so`). An ocean-only
reference serves the ocean of the rembo climate; `climate = "snapshots"`
needs the atmosphere. Checked at start-up, with an error naming the group:

- every snapshot of an axis of a `linear`/`ratio` field binds that field, and
  no snapshot binds a field the reference does not;
- a snapshot with `tas` or `pr` binds `zs` (the sea-level reduction), except
  the `hi` of a `lo = ""` axis (below); `pr` with `zs` needs `tas`; `so` needs
  `to`; the reference binds `tas` and `pr` together;
- the number of groups: `tas` one (`monthly = True`) or two (annual mean and one
  season); `pr` one or two summed (monthly, e.g. snowfall + rainfall) or one
  (annual); every other field one.

**Seasons.** A snapshot with two `tas` groups (annual mean `ta_ann` and one
season `ta_seas`) builds the monthly cycle as a cosine with its extremes in
January and July, `c(m) = cos(2π(m-1)/12)`: with `tas_season = "djf"`,
`tas(m) = ta_ann + (ta_seas - ta_ann) c(m)` (the season in January); with
`"jja"`, `tas(m) = ta_ann - (ta_seas - ta_ann) c(m)` (in July). The default is
the summer of the hemisphere (`jja` north, `djf` south), so the second field of
CLIMBER-3α (JJA only) serves the south as its winter, and PMIP3 offers both
(`t2m_jja`, `t2m_djf`).

**Ocean.** `to` and `so` are depth profiles, interpolated to the model depth
axis (23 levels, 0-3000 m). A monthly climatology (12 slices, e.g. ORAS5) is
averaged to the annual mean; an annual field (one slice: no time axis, e.g.
the CLIMBER ocean products, or one record) is taken as is. Without `so`, the
salinity is 0.

Each axis of a field adds to the reference snapshot, with weight `a_k` (the
index value):

| blend | now |
|---|---|
| `linear` | `ref + sum_k a_k (hi_k - lo_k)` |
| `ratio` | `ref * prod_k (1 + a_k (hi_k/lo_k - 1))` |
| `const` | `ref` |

An axis with `lo = ""` adds its `hi` snapshot as the anomaly itself:
`ref + a_k hi_k` (`linear`) or `ref (1 + a_k (hi_k - 1))` (`ratio`), e.g. a
reconstructed temperature anomaly or precipitation ratio (Badgeley et al.,
2020) without a zero snapshot. Such a snapshot binds no `zs`: its fields are
taken as given, so they must be **sea-level** anomalies or ratios (there is
no elevation correction of their own). If it binds `zs`, it is reduced to sea
level like any other snapshot.

```fortran
&snap_axis_b20                ! yelmox_GreenlandPaleo_deglac.nml
    lo    = ""
    hi    = "b20"             ! snapshot binding tas (anomaly) and pr (ratio), no zs
    index = "one"             ! a const index, f = 1
/
```

Two axes give the orbital + millennial decomposition of Banderas et al. (2018), e.g. `at` on
`picontrol` -> `lgm` and `bt` on `lgm_weak` -> `lgm_strong`. Several axes or
fields may share one index:

```fortran
&indices
    names = "at" "bt"
/

&idx_at
    method      = "series"
    series_file = "input/alpha_orbital_r2020.dat"
    norm_t0     = 0.0             ! [yr] index 0 here
    norm_t1     = -21e3           ! [yr] index 1 here (or a window, e.g. -23e3, -19e3)
/

&idx_bt
    method      = "series"
    series_file = "input/beta_ngrip_millennial.dat"
/

&snap_axis_at
    lo    = "picontrol"
    hi    = "lgm"
    index = "at"
/

&snap_axis_bt
    lo    = "lgm_weak"
    hi    = "lgm_strong"
    index = "bt"
/

&snap_field_tas
    axes  = "at" "bt"             ! blend = "linear" (default)
/
```

Turning a proxy series into an index is done by tsgen (`norm_t0`/`norm_t1`, see
the fesm-utils tsgen docs); a steady state is a `const` index (`df_dt_max = 0`,
value `f_min`), independent of the model time.

**The ocean anomaly from the atmosphere** is the climate modifier
`[ocn_fraction]` (see [yelmox](yelmox.md)): the field `to` without axes (or
`blend = "const"`) stays at the reference, and the modifier adds `f_to` times
the mean sea-level air temperature anomaly of the same update, without lag and
without the lapse-rate cooling of the ice sheet. `blend = "fraction"` stops the
run.

The backend holds no state of its own: the indices restart with the driver
(`idx_<name>_restart.nc`). The climate output carries the position on each
axis, `axis_<name>`.

The transient forcing is added on top in every configuration.
