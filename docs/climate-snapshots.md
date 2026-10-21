---
title: "Snapshots"
aliases:
  - climate-snap.html
---

The snapshot backend `snapshots` builds the climate from a reference
**snapshot** (e.g. present day) plus the **terms** of `[forcing]`: each a
scalar **series** on the forcing time (e.g. a glacial index) times a pattern,
the difference of two snapshots (e.g. LGM minus pre-industrial) or a uniform
anomaly (see [yelmox](yelmox.md#forcing)),
on the fesm-utils `varslice` reader and `tsgen` series. It supplies the
atmosphere (monthly `tas`, `pr`, lapse-rate corrected to the current surface)
and the ocean as depth profiles (`to_ann`, `so_ann`), which the marine shelf
interpolates to the shelf base. It runs on `grid_clim`, reads `[snapshots]`
and the `[snapshot_<name>]` groups, and also supplies the ocean of the rembo
climate (an ocean-only reference snapshot).

- **Program:** `yelmox` (`make yelmox`), backend in `libs/snapshots.f90`, behind
  `libs/yelmox_climate.f90`.
- **Select:** `[comps] climate = "snapshots"`.
- **Output:** with `[output] write_clim`, `snap.nc`, the climate as applied
  (after the transforms): the annual air temperature and its summer mean
  (`t2m_ann`, `t2m_sum`), the annual precipitation, the ocean profiles
  (`to_ann`, `so_ann` on `depth`), the value of each term's series
  (`term_<name>`), the surface elevation, the monthly atmosphere (`tas`,
  `pr`), `tsl_ann` and the diagnostics of the active transforms.

The transforms `[taper]` and `[bias]` (see [yelmox](yelmox.md)) apply after
the terms, in this order.

## Configuration

Example: `yelmox/yelmox_Greenland.nml`, present day plus the PMIP3 LGM minus
piControl anomaly weighted by `alpha` (here the orbital record r2020), the ocean plus a quarter of
the mean air temperature anomaly (the derived series `dTo`):

```fortran
&forcing
    reference = "pd"                  ! [snapshot_pd], in the par file
    series    = "alpha" "dTo"         ! [series_alpha] (par file: alpha_orbital_r2020), [series_dTo] (library)
    terms     = "lgm" "ocean"         ! [term_lgm]: picontrol -> lgm on tas, pr; [term_ocean]: dTo on to
/

&snapshots
    lapse     = 0.0080 0.0065         ! [K m-1] annual, summer (default)
    f_p       = 0.06                  ! [K-1] precipitation scaling (default 0.05)
/
```

The snapshots of the run are the reference and the `snap0`/`snap1` of its
terms (no list): `pd`, `picontrol` and `lgm` here. Every key of `[snapshots]`
and `[snapshot_<name>]` has a default in `input/yelmox_defaults.nml` (see
[Defaults](#defaults)): `combine = "anomaly"` (snapshots as anomalies from the
reference), `fields = "tas" "pr" "to" "so" "zs"` (the fields a snapshot may
bind), `grid_global = "GLOBAL-1DEG"`, ... Two kinds of group:

- **Rules**, `[snapshots]`: the lapse rates `lapse` (annual, summer) and the
  precipitation sensitivity `f_p` (β_p), used in the sea-level reduction of
  the snapshots, the uniform terms and the move to the model surface.
- **Snapshots**, `[snapshot_<name>]`: `exp` (the experiment, `{exp}` in the
  database paths), `monthly` (default `True`), `clim_sigma`/`ocn_sigma` ([km]
  Gaussian smoothing of the atmosphere and the ocean at loading, the land mask
  following the smoothed elevation; default 0), `tas_season` and, per field,
  the database group(s) holding it (see [Database](#database)).

### Database

The snapshots bind varslice groups (see [Database namelists](database-namelists.md))
by name, from the par file itself if it has the group (searched implicitly
first, so it never has to name itself and its groups win), else from the
files listed in `[forcing] var_defs`, searched in order (first file that has
the group). The snapshot groups `[snapshot_<name>]`, the series and the terms
are looked up the same way. The shared database `input/clim_db.nml` holds one
group per product field; a par file can add run-specific groups:

```fortran
&forcing
    var_defs = "input/clim_db.nml" "input/forcing_db.nml"   ! the default
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

- `{grid_global}`: the global grid of the paleoclimate products, `[snapshots]
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
| `Paleoclimate-CLIMBERX-Willeit2025-{exp}`, `-ocean` | piControl, interstadial, stadial, heinrich | `cx_{tas_mon,pr_mon,tas_ann,pr_ann,zs}`, `cx_{tas_mon_anom,pr_mon_anom}` (anomalies to piControl, for terms without `snap0`; `pr_mon_anom` is a difference in mm d-1, not a ratio), `cx_{to,so}` |

The Greenland reconstructions are regional products on GRL-16KM, remapped
online (`grid_src = "GRL-16KM"`) and filled from the nearest cell outside
their box (`fill = "nearest"`, e.g. on GRL-PAL):

| Product (FesmData) | Time axis | Groups |
|---|---|---|
| `Paleoclimate-Buizert2018` | monthly, every 10 yr, -22005 to 55 yr (22005 BP to 2005 CE) | `b18_{tas,pr}` (time-varying), `b18_zs` (static) |
| `Paleoclimate-Buizert2018-1981-2010` | none (its 1981-2010 baseline) | `b18_1981_2010_{tas,pr,zs}` |
| `Paleoclimate-Badgeley2020` | annual, every 50 yr, 20 ka to 0 | `b20_tas_anom` [K], `b20_pr_frac`, `b20_pr_frac_low`, `b20_pr_frac_high` (anomaly and ratios relative to 1850-2000 CE, for terms without `snap0`) |

The database also holds the groups of the transforms (`[bias]`:
`bias_mar_era5_1981_2010_{dtas,fpr}`; `[taper]`: `taper_b18_b20`) and the ice
histories of `[ice_bnd]` (`paleomist_H_ice`, `glac1d_H_ice`).

**Time-varying groups.** A group with a time axis (`time = 1 x0 x1 dx nsub`)
is sliced at the forcing time `time_rel` at the start and at every climate
update: linear interpolation between the two bracketing records, the first or
last record outside the axis (varslice `extrap`). Its snapshot is then reduced
to sea level again (and the reference reset, if it is the reference). A group
without a time axis (`time = 0 ...`) is read once at the start.

```fortran
&snapshot_pd
    tas = "mar314_1961_1990_tas"
    pr  = "mar314_1961_1990_sf" "mar314_1961_1990_rf"
    zs  = "mar314_1961_1990_zs"
    to  = "oras5_1981_2010_to"
    so  = "oras5_1981_2010_so"
/

&snapshot_lgm
    exp        = "lgm"
    monthly    = False
    clim_sigma = 50.0
    tas = "pmip3_mean_tas_ann" "pmip3_mean_tas_jja"
    pr  = "pmip3_mean_pr_ann"
    zs  = "pmip3_mean_zs"
/
```

### Defaults

`input/yelmox_defaults.nml` declares every key of `[snapshots]` and
`[snapshot_<name>]`, of the forcing (`[forcing]`, the series, the terms) and
of the transforms, with its default. The named groups fall back to a template
group: `[snapshot_<name>]` to `[snapshot]`, `[snapshots<sfx>]` to
`[snapshots]` with the domain suffix of a multi-domain driver
(`[snapshots_north]`). A par file lists only what differs; a key the defaults
file does not declare stops the run (`nml_validate`), and
`scripts/check_par_nml.py` reports it offline.

### Fields and composition

**Fields are optional per snapshot.** A snapshot binds the fields it names
(a key absent or `""` binds nothing). The backend builds what the reference
snapshot (`[forcing] reference`) binds: the atmosphere (`tas` and `pr`, with
the elevation `zs`) and/or the ocean (`to`, optionally `so`). An ocean-only
reference serves the ocean of the rembo climate; `climate = "snapshots"`
needs the atmosphere. Checked at start-up, with an error naming the group:

- every snapshot of a term binds the fields the term lists, and no snapshot
  binds a field the reference does not;
- a snapshot with `tas` or `pr` binds `zs` (the sea-level reduction), except
  the `snap1` of a term without `snap0` (below); `pr` with `zs` needs `tas`;
  `so` needs `to`; the reference binds `tas` and `pr` together;
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

**Composition.** Each term with the `snapshots` pattern adds to the reference
with the weight `a_k`, the value of its series; the rule is fixed by field:

| field | now |
|---|---|
| `tas`, `to`, `so` | `ref + sum_k a_k (snap1_k - snap0_k)` |
| `pr` | `ref * prod_k (1 + a_k (snap1_k/snap0_k - 1))` |

A term with `snap0 = ""` adds its `snap1` snapshot as the anomaly itself:
`ref + a_k snap1_k` (`tas`) or `ref (1 + a_k (snap1_k - 1))` (`pr`), e.g. a
reconstructed temperature anomaly or precipitation ratio (Badgeley et al.,
2020) without a zero snapshot. Such a snapshot binds no `zs`: its fields are
taken as given, so they must be **sea-level** anomalies or ratios (there is
no elevation correction of their own). If it binds `zs`, it is reduced to sea
level like any other snapshot.

```fortran
&term_b20                     ! input/forcing_db.nml
    snap0  = ""
    snap1  = "b20"            ! snapshot binding tas (anomaly) and pr (ratio), no zs
    fields = "pr"             ! series = "": the constant 1
/
```

Two terms give the orbital + millennial decomposition of Banderas et al.
(2018), e.g. `lgm` on `picontrol` -> `lgm` and `millennial` on
`lgm_weak` -> `lgm_strong`, each on its own series; several terms may share
one series:

```fortran
&forcing
    reference = "pd"
    series    = "alpha" "beta"
    terms     = "lgm" "millennial"
/

&series_alpha
    method      = "series"
    series_file = "input/alpha_orbital_r2020.dat"
    norm_t0     = 0.0             ! [yr] series 0 here
    norm_t1     = -21e3           ! [yr] series 1 here (or a window, e.g. -23e3 -19e3: the mean)
/

&series_beta
    method      = "series"
    series_file = "input/beta_ngrip_millennial.dat"
/

&term_lgm
    series = "alpha"
    snap0  = "picontrol"
    snap1  = "lgm"
    fields = "tas" "pr"
/

&term_millennial
    series = "beta"
    snap0  = "lgm_weak"
    snap1  = "lgm_strong"
    fields = "tas" "pr"
/
```

Turning a proxy series into an index is done by tsgen (`norm_t0`/`norm_t1`,
see the fesm-utils tsgen docs); a steady state is a `const` series
(`df_dt_max = 0`, value `f_min`), independent of the model time.

**Uniform terms** (`pattern = "uniform"`) add `f_ta a` to `tas` (and multiply
`pr` by exp(β_p f_ta a)), `f_to a` to `to` and `f_so a` to `so` on the fields
they list. **The ocean anomaly from the atmosphere** is a derived series
(`method = "derived"`, `from = "air_anomaly"`: `factor` times the mean
sea-level air temperature anomaly of the same update over `region`, without
lag and without the lapse-rate cooling of the ice sheet) with a uniform term
on `to`, the library `series_dTo` + `term_ocean`; `to` without terms stays at
the reference. The ocean terms are applied after the atmosphere (a derived
series may drive ocean fields only).

The backend holds no state of its own: the series restart with the driver
(`series_<name>_restart.nc`). The climate output carries the value of each
term's series, `term_<name>`.

The transient forcing is added on top in every configuration.
