---
title: "snapclim and snapesm"
---

The two snapshot backends build the climate from a reference state plus climate
**snapshots** (e.g. present day, pre-industrial, LGM), blended in time by scalar
**indices** (e.g. a glacial index). Both supply the atmosphere (monthly `tas`,
`pr`, lapse-rate corrected to the current surface) and the ocean as depth
profiles (`to_ann`, `so_ann`), which the marine shelf interpolates to the shelf
base. Both run on `grid_clim` and read the `[snap]` group. snapesm also
supplies the ocean of the rembo climate (an ocean-only reference snapshot).

- **Program:** `yelmox` (`make yelmox`), backends in `libs/snapclim.f90` and
  `libs/snapesm.f90`, behind `libs/yelmox_climate.f90`.
- **Select:** `[comps] climate = "snapclim"` or `"snapesm"`.
- **Output:** with `[output] write_clim`, `snap.nc`, the climate as applied
  (after the climate modifiers): the annual air temperature and its summer
  mean (`t2m_ann`, `t2m_sum`), the annual precipitation and the ocean profiles
  (`to_ann`, `so_ann` on `depth`); with snapesm also the position on each axis
  (`axis_<name>`), the surface elevation, the monthly atmosphere (`tas`, `pr`)
  and `tsl_ann`.

## snapclim

The default backend of every regional config (`yelmox/yelmox_<domain>.nml`).
`[snap]` sets the forcing method for the atmosphere (`atm_type`) and the ocean
(`ocn_type`), the index files (`fname_at`, `fname_ao`, ...), the lapse rates and
the scaling parameters; up to four snapshots are read from `[snap_clim0]` to
`[snap_clim3]`, plus `[snap_hybrid]` and `[snap_recon]` for those methods.

| Method | Climate |
|---|---|
| `const` | the reference snapshot (`clim0`) |
| `anom` | reference + a homogeneous anomaly (`dTa`/`dTo`/`dSo`) from its index series |
| `snap_1ind`, `snap_1ind_new`, `snap_1ind_abs`, `snap_1ind_miocene` | two snapshots blended by one index |
| `snap_2ind`, `snap_2ind_abs` | snapshots blended by two indices |
| `hybrid` | reference + monthly temperature anomalies from a prescribed series (`[snap_hybrid]`), no snapshots |
| `recon` | reconstructed snapshots, interpolated in time (`[snap_recon]`) |
| `fraction` (ocean only) | an ocean anomaly of `f_to` × the mean atmospheric anomaly |

The climate modifiers `[anom]` and `[ocn_fraction]` (see [yelmox](yelmox.md))
apply after the backend, to any method.

A snapshot group reads the atmosphere (`clim_path`, `clim_names`, `clim_monthly`,
`clim_stdev_path`/`clim_stdev_name`) and the ocean (`ocn_path`, `ocn_names`,
`ocn_monthly`) of one snapshot, from the [FesmData](https://github.com/fesmc/FesmData)
v2 products (`fesmdata/`), e.g. the Greenland reference and LGM snapshot:

```fortran
&snap_clim0
    clim_path       = "fesmdata/{domain}/{grid_name}/{grid_name}_SMB-MARv3.14-1961-1990.nc"
    clim_names      = "z_srf" "t2m" "sf" "rf"
    clim_monthly    = True
    ...
    ocn_path        = "fesmdata/{domain}/{grid_name}/{grid_name}_Ocean-ORAS5-1981-2010.nc"
    ocn_names       = "depth" "mask_ocn" "to" "so"
    ocn_monthly     = True
/
&snap_clim2
    clim_path       = "fesmdata/{domain}/{grid_name}/{grid_name}_Paleoclimate-PMIP3-lgm-mean.nc"
    clim_names      = "z_srf" "t2m_ann" "t2m_jja" "pr_ann"
    clim_monthly    = False
    clim_time       = -21000
    clim_stdev_path = "fesmdata/{domain}/{grid_name}/{grid_name}_Paleoclimate-PMIP3-lgm-stdev.nc"
    clim_stdev_name = "pr_ann"
    clim_sigma      = 50.0      ! [km] smoothing of the atmosphere
    ocn_path        = "fesmdata/{domain}/{grid_name}/{grid_name}_Paleoclimate-PMIP3-lgm-mean.nc"
    ocn_names       = "none" "none" "sst_ann"
    ocn_monthly     = False
    ocn_time        = -21000
    ocn_sigma       = 50.0      ! [km] smoothing of the ocean
/
```

- **Units** come from the `units` attribute of each variable: temperatures in K
  or °C, precipitation in any water-flux unit (kg m⁻² or mm w.e. per s, d or yr);
  an unknown unit stops the run.
- **Smoothing.** `clim_sigma` and `ocn_sigma` (km, optional, default 0) smooth the
  snapshot's atmosphere and ocean with a Gaussian filter when it is loaded; the
  land mask follows the smoothed elevation. The PMIP3 snapshots use 50 km (the v1
  files were smoothed in preprocessing, `PMIP3_sig50km`).
- **Ocean gaps.** Cells without ocean data (land, or beyond the coverage of the
  product) are filled level by level from the nearest valid cell, and a level
  without data from the level above, so the ocean is defined wherever the ice
  sheet meets it. The ocean mask (`ocn_names(2)`) is read from the file, or
  taken from where the data exist with `"none"`. A surface ocean
  (`ocn_names(1) = "none"`, e.g. PMIP3 `sst_ann`) has salinity 0.

## snapesm

A rewrite of snapclim on the fesm-utils `varslice` reader and `tsgen` series:
the snapclim methods become configurations of one blend model, over any number of
snapshots, fields and axes. Example: `yelmox/yelmox_Greenland_snapesm.nml`,
which reproduces the atmosphere of snapclim `snap_1ind_new`
(`yelmox_Greenland_snapclim.nml`); its ocean anomaly comes from `[ocn_fraction]`.

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
  default `linear`) and the `axes`, one or more names (default `""` = none),
  and `kind` (`atm_temp`, `atm_precip`, `ocn`, `ocn_salt`, `elev`; default `""`,
  the kind of the field name). A field without axes follows the reference: its
  group can be left out.
- **Snapshots**, `[snap_<name>]`: `exp` (the experiment, `{exp}` in the
  database paths), `monthly` (default `True`), `clim_sigma`/`ocn_sigma` (as in
  snapclim, default 0), `tas_season` and, per field, the database group(s)
  holding it (see [Database](#database)).
- **Axes**, `[snap_axis_<name>]`: two snapshots `lo` -> `hi` and the `index`
  giving the position on the axis: 0 at `lo`, 1 at `hi`. The index is a
  driver-owned `[idx_<name>]` (a `tsgen` generator of any method: `const`,
  ramps, `series`, ...; see [yelmox](yelmox.md)), listed in `[indices]`.

### Database

The snapshots bind varslice groups (see [Database namelists](database-namelists.md))
by name, from the files listed in `[snap] var_defs`, searched in order (first
file that has the group), then from the par file itself (searched implicitly
last, so it never has to name itself). The shared database `input/clim_db.nml`
holds one group per product field; a par file can add run-specific groups:

```fortran
&snap
    var_defs = "input/clim_db.nml"   ! the default
/
&my_tas                              ! in the par file: found after clim_db.nml
    filename = "ice_data/..."
    ...
/
```

A group missing from every file (and the par file) stops the run, naming the
group and the files.
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
| `Paleoclimate-CLIMBERX-Willeit2025-{exp}`, `-ocean` | piControl, interstadial, stadial, heinrich | `cx_{tas_mon,pr_mon,tas_ann,pr_ann,zs}`, `cx_{tas_mon_anom,pr_mon_anom}` (anomalies to piControl, for `lo = ""` axes), `cx_{to,so}` |

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
reference serves the ocean of the rembo climate; `climate = "snapesm"` needs
the atmosphere. Checked at start-up, with an error naming the group:

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
salinity is 0 (as in snapclim).

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
&snap_axis_recon
    lo    = ""
    hi    = "badgeley"      ! snapshot binding tas (anomaly) and pr (ratio), no zs
    index = "one"
/
```

One axis gives snapclim `snap_1ind`; two axes give `snap_2ind`, the
orbital + millennial decomposition of Banderas et al. (2018), e.g. `at` on
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

**The ocean anomaly from the atmosphere** (snapclim `fraction`) is the climate
modifier `[ocn_fraction]` (see [yelmox](yelmox.md)): the field `to` without
axes (or `blend = "const"`) stays at the reference, and the modifier adds `f_to` times
the mean sea-level air temperature anomaly of the same update. Unlike
snapclim's `fraction` (and snapesm's before), it has no lag of one update and
leaves out the lapse-rate cooling of the ice sheet. `blend = "fraction"` stops
the run.

snapesm holds no state of its own: the indices restart with the driver
(`idx_<name>_restart.nc`). The climate output carries the position on each
axis, `axis_<name>`.

The transient forcing is added on top in every configuration.
