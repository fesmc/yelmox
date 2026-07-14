---
title: "Database namelists (varslice)"
---

Most forcing in YelmoX is now loaded through the shared
[`varslice`](https://github.com/fesmc/fesm-utils) reader. A **database namelist**
is a small, static description of a dataset on disk: where the file lives, which
variable to read, how to convert its units, and how its records map onto a time
axis. The consuming modules — `varslice` itself, `esm_forcing`, `snapesm`,
`ismip6` — all read the *same* group layout, so once you understand one group you
understand them all.

The point of these files is to keep the *provenance* of every forcing field
visible and in one place. Nothing in the namelist runs code; each group is a
declarative record that says "this NetCDF file, this variable, these units, this
period" — everything the reader needs to ingest the dataset and hand back a
time-indexed field.

See `input/esm/esm_ant_ismip7.nml` for a full, real-world example.

## Anatomy of a group

A single dataset is described by one namelist group. All five keys are
**mandatory** (the reader pre-sets defaults, but the convention is to always
write them out so the provenance stays in the file):

```fortran
&gcm_ts_ref
    filename    = "ice_data/{domain}/{grid_name}/RACMO2.3/{grid_name}_ERA5-3H_RACMO2.3p2_1979-2022_monthly.nc"
    name        = "t2m"
    units       = "K" "K"
    scaling     = 1.0 0.0
    time        = 1, 1979.0, 2022.0, 1.0, 12.0
/
```

| Key        | Type            | Meaning |
|------------|-----------------|---------|
| `filename` | string          | Path to the NetCDF file (relative to the run dir). May contain placeholders and wildcards — see [Paths](#paths-placeholders-and-wildcards). |
| `name`     | string          | Name of the variable to read inside the NetCDF file. |
| `units`    | 2 strings       | `"<units_in>" "<units_out>"` — documentation labels only. |
| `scaling`  | 2 reals         | `<unit_scale> <unit_offset>` — the actual linear unit conversion. |
| `time`     | 5 reals         | `<active> <x0> <x1> <dx> <sub>` — describes the time axis. See [The `time` array](#the-time-array). |

The group **name** (`&gcm_ts_ref`, `&gcm_to_proj`, …) is chosen by the consuming
module — it decides which physical role a group plays (reference vs. projection,
atmosphere vs. ocean, etc.). The five keys inside are always the same.

### `units` vs. `scaling`

The two `units` strings are **labels for humans** — they do not trigger any
conversion. The real conversion is the linear transform in `scaling`, applied to
every value as it is read:

```
value_out = value_in * unit_scale + unit_offset
```

`units` documents the intent so a reader can check that `scaling` actually
performs it. For example, precipitation stored in `kg m-2 s-1` and wanted in
`mm yr-1`:

```fortran
units   = "kg m-2 s-1" "mm yr-1"
scaling = 31556926.25 0.0     ! CMIP standard seconds-per-year
```

or ocean temperature stored in `°C` and wanted in `K`:

```fortran
units   = "degrees C" "K"
scaling = 1.0 273.15
```

Use `scaling = 1.0 0.0` when no conversion is needed.

## The `time` array

The `time` key is the only subtle one. It packs five numbers:

```
time = active, x0, x1, dx, sub
```

| Element  | Name     | Meaning |
|----------|----------|---------|
| `active` | with_time| `1` = the field is **time-varying** (read record-by-record). `0` = the field is **static** (one slice); `x0..x1` then merely *document* the period the data represents. |
| `x0`     | start    | First time value on the axis (calendar year, in ESM convention). |
| `x1`     | end      | Last time value on the axis. |
| `dx`     | step     | Spacing between successive **major** time units (usually `1.0` year). |
| `sub`    | sub-count| Number of **sub-annual** records per major unit: `1` = annual, `12` = monthly. Any value `> 1` switches on sub-annual handling. |

How the reader interprets these:

- **Static field** (`active = 0`). `x0..x1` are preserved as-is and describe the
  period the (single) slice represents. Typical for a fixed surface elevation or
  an ocean climatology:

  ```fortran
  time = 0, 2000, 2000, 0, 1     ! one slice, "representing year 2000"
  time = 0, 1990, 2014, 0, 1     ! one slice, "a 1990–2014 climatology"
  ```

- **Time-varying, annual** (`active = 1`, `sub = 1`). The axis is built as
  `x0, x0+dx, …, x1`:

  ```fortran
  time = 1, 1850, 2014, 1, 1     ! 165 annual records, 1850…2014
  ```

- **Time-varying, sub-annual** (`active = 1`, `sub > 1`). A sub-annual axis is
  built with spacing `1/sub`, centred within each major unit (e.g. mid-month for
  `sub = 12`), spanning `x0` through `x1`:

  ```fortran
  time = 1, 1979.0, 2022.0, 1.0, 12.0   ! monthly records, Jan 1979 … Dec 2022
  ```

- **Static field, sub-annual** (`active = 0`, `sub > 1`). `active` and `sub` are
  independent switches: `active` decides whether the field advances along a
  multi-year calendar axis, `sub` decides whether there is a within-year
  (monthly) cycle. Combining them gives a **monthly climatology** — a fixed
  seasonal cycle, selected by month-of-year, that does *not* advance year to year.
  `x0..x1` document the period the climatology averages over; `dx = 0`:

  ```fortran
  time = 0, 1961, 1990, 0, 12     ! 12-month climatology, averaged over 1961–1990
  ```

- **Collapse to a single slice.** For a *time-varying* field, setting `dx = 0`
  collapses `x1` to `x0`, i.e. read exactly one record at `x0`. (For a static
  field, `dx = 0` is normal and `x0..x1` are left untouched.)

The generated axis length must match the length of the NetCDF file's time
dimension, or the reader stops with an error — a useful check that `x0`, `x1`,
`dx`, `sub` describe the file you actually pointed at.

## Paths: placeholders and wildcards

`filename` is expanded before the file is opened.

**Placeholders** in `{braces}` are substituted by the caller:

- `{domain}` and `{grid_name}` are built in — every module supplies them (e.g.
  `Antarctica`, `ANT-16KM`).
- Additional keys are passed by the consuming module. `esm_forcing`, for
  instance, supplies `{gcm}` and `{experiment}`, so an ISMIP7 projection path
  can read:

  ```fortran
  filename = "ice_data/ISMIP7/{domain}/{grid_name}/{gcm}/{experiment}/ocean/thetao/v3/{grid_name}_{gcm}_{experiment}_v3_2015-2299_thetao.nc"
  ```

  and the same namelist serves every GCM/experiment combination — you select the
  concrete values in the driver's `[esm]` configuration, not here.

**Wildcards** (shell globs, `*` and `?`) are also allowed. The reader expands the
pattern with `ls -1 … | sort`, so a dataset split across many files is matched
and concatenated in alphabetical order:

```fortran
filename = "ice_data/{domain}/{grid_name}/{gcm}/historical/tas/*.nc"
```

Sort the split so alphabetical order equals chronological order (zero-pad years).

## Worked examples

Monthly, time-varying reference climate (RACMO2.3 on the ice grid):

```fortran
&gcm_ts_ref
    filename    = "ice_data/{domain}/{grid_name}/RACMO2.3/{grid_name}_ERA5-3H_RACMO2.3p2_1979-2022_monthly.nc"
    name        = "t2m"
    units       = "K" "K"
    scaling     = 1.0 0.0
    time        = 1, 1979.0, 2022.0, 1.0, 12.0
/
```

Annual ESM projection ocean temperature, `°C → K`, GCM/experiment templated:

```fortran
&gcm_to_proj
    filename    = "ice_data/ISMIP7/{domain}/{grid_name}/{gcm}/{experiment}/ocean/thetao/v3/{grid_name}_{gcm}_{experiment}_v3_2015-2299_thetao.nc"
    name        = "thetao"
    units       = "degrees C" "K"
    scaling     = 1.0 273.15
    time        = 1, 2015, 2299, 1, 1
/
```

Static surface elevation (single slice, no unit change):

```fortran
&gcm_zs_ref
    filename    = "ice_data/{domain}/{grid_name}/RACMO2.3/{grid_name}_ERA5-3H_RACMO2.3p2_1979-2022_monthly.nc"
    name        = "z_srf"
    units       = "m" "m"
    scaling     = 1.0 0.0
    time        = 0, 2000, 2000, 0, 1
/
```

Static ocean climatology (single slice representing a 1990–2014 average):

```fortran
&gcm_to_ref
    filename    = "ice_data/ISMIP7/{domain}/{grid_name}/meltMIP/{grid_name}_OI_Climatology_thetao_extrap.nc"
    name        = "thetao"
    units       = "degrees C" "K"
    scaling     = 1.0 273.15
    time        = 0, 1990, 2014, 0, 1
/
```

## How the groups are consumed

A module loads a group into a `varslice` object once at initialization
(`varslice_init_nml(vs, filename, group, domain, grid_name, subs=…)`) and then
asks for a field at a given time during the run (`varslice_update(vs, time,
method=…)`). The database namelist supplies everything the object needs; the
`method` (e.g. `exact`, `interp`, `extrap`) is chosen by the caller at read time,
not in the namelist.

- **`esm_forcing`** groups the atmosphere/ocean database groups into reference,
  variability, ESM-reference, historical, and projection periods
  (`&gcm_ts_ref`, `&gcm_to_proj`, …) and supplies `{gcm}`/`{experiment}`. See
  [yelmox_esm](flavor-esm.md).
- **`snapesm`** uses the same reader for its snapshot climate database.
- **`ismip6`** uses it for the ISMIP6 atmosphere/ocean forcing files.

Because they share the reader, the layout documented here is identical across all
of them — a group is a group regardless of which module owns it.

## Quick reference

- All five keys are mandatory; always write them for provenance.
- `units` = two documentation labels; `scaling` = `[scale, offset]` does the real
  `out = in*scale + offset` conversion.
- `time = [active, x0, x1, dx, sub]`:
  `active` 1/0 = time-varying/static · `sub` 1/12 = annual/monthly ·
  `dx = 0` collapses a time-varying field to a single slice at `x0`.
- Paths expand `{domain}`/`{grid_name}` plus module-supplied keys like
  `{gcm}`/`{experiment}`, and accept `*`/`?` globs (sorted alphabetically).
- The generated time axis length must match the file's time dimension.
