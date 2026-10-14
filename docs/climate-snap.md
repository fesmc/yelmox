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
  position on the index manifold), optionally `clim_sigma`/`ocn_sigma` (as in
  snapclim) and, per field, the `var_defs` group(s) holding it (one = monthly;
  two = annual + summer, from which the seasonal cycle is built). See
  [Database namelists](database-namelists.md).
- **Indices**, `[snap_idx_<name>]`: a `tsgen` series (`method`, `series_file`,
  `sigma`).

Restart bundles hold the state of each index, `snapesm_idx_<name>_restart.nc`,
with the configuration (snapshots, fields, reference) as attributes.

The transient forcing is added on top in every configuration.
