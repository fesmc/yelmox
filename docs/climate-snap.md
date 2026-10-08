---
title: "snapclim and snapesm"
---

The two snapshot backends build the climate from a reference state plus climate
**snapshots** (e.g. present day, pre-industrial, LGM), blended in time by scalar
**indices** (e.g. a glacial index). Both supply the atmosphere (monthly `tas`,
`pr`, lapse-rate corrected to the current surface) and the ocean as depth
profiles (`to_ann`, `so_ann`), which the marine shelf interpolates to the shelf
base. Both run on `grid_clim` and read the `[snap]` group.

- **Program:** `yelmox` (`make yelmox`), backends in `libs/snapclim.f90` and
  `libs/snapesm.f90`, behind `libs/yelmox_climate.f90`.
- **Select:** `[comps] climate = "snapclim"` or `"snapesm"`.
- **Output:** with `[output] write_clim`, `snap.nc`: the annual air temperature
  and precipitation; with snapesm also its driving indices (`idx_<name>`), the
  surface elevation, the monthly atmosphere (`tas`, `tsl`, `pr`), `ta_sum`,
  `tsl_ann` and the ocean profiles (`to_ann`, `so_ann` on `depth`).

## snapclim

The default backend of every regional config (`yelmox/yelmox_<domain>.nml`).
`[snap]` sets the forcing method for the atmosphere (`atm_type`) and the ocean
(`ocn_type`), the index files (`fname_at`, `fname_ao`, ...), the lapse rates and
the scaling parameters; up to four snapshots are read from `[snap_clim0]` to
`[snap_clim3]`, plus `[snap_hybrid]` and `[snap_recon]` for those methods.

| Method | Climate |
|---|---|
| `const` | the reference snapshot (`clim0`) |
| `anom` | reference + a homogeneous anomaly (`dTa`/`dTo`/`dSo`): the transient forcing's, else from the index |
| `snap_1ind`, `snap_1ind_new`, `snap_1ind_abs`, `snap_1ind_miocene` | two snapshots blended by one index |
| `snap_2ind`, `snap_2ind_abs` | snapshots blended by two indices |
| `hybrid` | reference + monthly temperature anomalies from a prescribed series (`[snap_hybrid]`), no snapshots |
| `recon` | reconstructed snapshots, interpolated in time (`[snap_recon]`) |
| `fraction` (ocean only) | an ocean anomaly of `f_to` × the mean atmospheric anomaly |

The transient forcing (`[tsforcing]`/`[tsgen]`, see [yelmox](yelmox.md))
reaches snapclim only in the `anom` methods.

## snapesm

A rewrite of snapclim on the fesm-utils `varslice` reader and `tsgen` series:
the snapclim methods become configurations of one blend model, over any number of
snapshots, fields and indices. Example: `yelmox/yelmox_Greenland_snapesm.nml`,
which reproduces snapclim `snap_1ind_new` / `fraction`
(`yelmox_Greenland_snapclim.nml`).

```fortran
&snap
    var_defs  = "input/greenland_clim.nml"   ! database of the input files
    combine   = "anomaly"                    ! snapshots as anomalies from the reference
    manifold  = 1                            ! number of indices spanning the blend (0/1/2)
    ref_name  = "pd"                         ! the reference snapshot
    snapshots = "pd" "picontrol" "lgm"
    fields    = "tas" "pr" "to" "so" "zs"
    indices   = "at" "ap"
    ! lapse, f_p, f_p_ne, f_stdev, f_to, f_hol, dTa_const, dTo_const, dSo_const
    ! as in snapclim
/
```

Three kinds of group follow:

- **Fields**, `[snap_field_<name>]`: `kind` (`atm_temp`, `atm_precip`, `ocn`,
  `ocn_salt`, `elev`), `blend` (`linear`, `ratio`, `fraction`, `const`) and the
  driving `index` (`""` = none).
- **Snapshots**, `[snap_<name>]`: `monthly`, `time`, `idx_coord` (the snapshot's
  position on the index manifold) and, per field, the `var_defs` group(s) holding
  it (one = monthly; two = annual + summer, from which the seasonal cycle is
  built). See [Database namelists](database-namelists.md).
- **Indices**, `[snap_idx_<name>]`: a `tsgen` series (`method`, `series_file`,
  `sigma`).

Restart bundles hold the state of each index, `snapesm_idx_<name>_restart.nc`,
with the configuration (snapshots, fields, reference) as attributes.

The transient forcing is added on top in every configuration.
