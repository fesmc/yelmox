---
title: "Greenland paleo"
---

A Greenland simulation of the last glacial period on the GreenlandPaleo domain
(GRL-PAL, 16 km), the Greenland domain extended 640 km to the west (Ellesmere,
Devon Island, Baffin Bay). In two steps, as the [Antarctic paleo](experiment-ant-paleo.md):
a 15-kyr present-day **spin-up** that optimizes the basal friction (`cb_ref`) and
the ocean thermal-forcing correction (`tf_corr`), then a **transient** from 130
kyr ago to the present, restarted from the spin-up, with the climate scaled by a
glacial index and the sea level following a reconstruction.

- **Program:** `yelmox` (`make yelmox`) with the snapshot backend
  (`[comps] climate = "snapshots"`, see [Snapshots](climate-snapshots.md)).
- **Config:** `yelmox/yelmox_GreenlandPaleo.nml` (base),
  `yelmox/yelmox_GreenlandPaleo_spinup.nml`, `yelmox/yelmox_GreenlandPaleo_lgp.nml`.
- **Scripts:** `scripts/grl-paleo/run_spinup.sh`, `scripts/grl-paleo/run_lgp.sh`.
- **Design:** `docs/_grl-paleo-plan.md` (setup A; the deglaciation after Tabone
  et al., 2024, setup B, is not set up yet).

## Running

Run from the yelmox root; the transient needs the spin-up's last restart bundle.

```bash
scripts/grl-paleo/run_spinup.sh
```

```bash
scripts/grl-paleo/run_lgp.sh
```

`run_spinup.sh [output_path]` writes to `output/grl-paleo/spinup` by default.
`run_lgp.sh [spinup_path] [output_path]` restarts from
`<spinup_path>/restart-15.000-kyr` (passed as an absolute path; the script stops if
the bundle is missing) and writes to `output/grl-paleo/lgp`. Both submit to the
`shared` queue (2 days, 10 GB, no OpenMP).

## Domain borders

The west border cuts Baffin Island and the Canadian Arctic Archipelago, where
Innuitian and Laurentide ice crosses it in glacial periods. The borders are
"infinite" (`[yelmo] experiment = "infinite"`, `mask_border = "dynamic"`): the
outermost cells take the thickness of their inner neighbours, so the border ice
thickens and thins with the interior. No ice enters from outside the domain: the
ice west of the border is what the local climate builds.

Prescribed border ice is the alternative (planned for the deglaciation, B): a
strip of fixed and relaxation cells (`[yelmo_masks] border_width = N M`, e.g.
`2 4`) whose thickness comes from an ice history (`[ice_bnd] active = True`,
PaleoMIST by default), with `experiment = "None"` and `mask_border = "auto"`.

## Reference climate

ERA5 1981--2010 (monthly) is the present-day reference on the whole domain; MAR
v3.14 covers only Greenland. Over Greenland, `[bias]` corrects ERA5 to MAR
1981--2010 with a monthly temperature difference `dtas` and precipitation factor
`fpr`, both taken at sea level with the reduction of the snapshot backend
(`[snap] lapse`, `f_p`), so they do not depend on elevation. Outside the MAR
coverage the correction is the nearest value of the bias smoothed over 50 km,
tapered to zero over 100 km (cosine). `fpr` is clamped to [0.2, 5].

The product is a hand-built file in ice_data,
`ice_data/GreenlandPaleo/GRL-PAL-16KM/GRL-PAL-16KM_BIAS-MARv3.14-ERA5-1981-2010.nc`
(`clim_db.nml` groups `bias_mar_era5_1981_2010_dtas`, `_fpr`), made by

```bash
julia --project=scripts/grl-paleo scripts/grl-paleo/make_bias_mar_era5.jl --fig
```

from the yelmox root (`fesmdata` and `ice_data` links in place; options `--res`,
`--l1`, `--l2`, `--f_p`, `--D`, `--sigma`, `--fpr_min`, `--fpr_max`). The lapse
rates and `f_p` are written as attributes: remake the file when `[snap] lapse` or
`f_p` change.

## Spin-up

- **Timeline:** `[ctrl]` 0 to 15 kyr, `dtt = 10` yr, present-day insolation
  (`[orbit] method = "time"`, `time = 0`), sea level from file at 0 ka.
- **Optimization** (`[sim] opt = True`): `cb_ref` (`L21`) and `tf_corr`
  (`L21-points`) over the whole 15 kyr; the ice thickness relaxes towards the
  observations until 6 kyr (`[sim] relax = True`, `[relax] time2`).
- **Climate:** present day (`[snap] snapshots = "pd"`) with the MAR -- ERA5 bias.
- **Output:** restart bundles at 0 and 15 kyr; 2D output every 1 kyr.

## Transient

- **Timeline:** `[ctrl]` -130 kyr to 0, `dtt = 10` yr, on relative time
  (`tstep_method = "rel"`). `[sim] relax = False`, `opt = False`: `cb_ref` and
  `tf_corr` come from the restart.
- **Climate:** snapshots `pd`, `picontrol`, `lgm` and the axis `a`
  (`[snap_axis_a]`): the present-day reference plus the PMIP3 LGM -- piControl
  anomaly, scaled by the glacial index `[idx_a]`
  (`input/alpha_combined_125kyr_interp.dat`, as is: 0 = present day, ~0.9 at
  21 ka, -0.2 in the last interglacial). The ocean anomaly is `f_to = 0.25` times
  the mean sea-level air temperature anomaly (`[ocn_fraction]`). Insolation varies
  in time (`[orbit] method = "clock"`).
- **Sea level:** `[barysealevel] method = "file"`, `input/sealevel_rohling_450kyr.dat`.
- **Ice age:** isochrones with elsa (`[ytrc] use_elsa`, `[elsa_grl_paleo]`): one
  layer every 1 kyr, coupled every 50 yr (for comparison with MacGregor et al., 2025).
- **Output:** 2D every 1 kyr, 1D every 100 yr, restart bundles every 25 kyr.

## Configuration common to both

| | |
|---|---|
| Grids | `GRL-PAL-16KM` for every component |
| Topography | BedMachine v6 (`GRL-PAL-16KM_TOPO-BedMachine-v6.nc`) |
| Reference climate (`snap_pd`) | ERA5 1981--2010 + MAR -- ERA5 bias over Greenland (`[bias]`); ocean ORAS5 1981--2010 |
| Snapshots (`snap_picontrol`, `snap_lgm`) | PMIP3 piControl and LGM ensemble means, global (`GLOBAL-1DEG`, remapped, 50 km smoothing; atmosphere only, the ocean follows `[ocn_fraction]`) |
| Surface mass balance | chion (`[comps] surface_method = "chion"`) |
| Shelf melt | `bmb_method = "lin"`, `tf_method = 1` |
| Dynamics | DIVA |
| Isostasy | FastIsostasy LV-ELRA, uniform mantle (the GRL-PAL Lloyd2024 file does not have the layout of `rheology_file`) |
| Borders | "infinite" (see above) |
