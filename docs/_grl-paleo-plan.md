# Greenland paleo: glacial cycles and a Tabone et al. (2024) deglaciation

Status: design summary (2026-10-17), not started. The leading `_` keeps
quarto from rendering it. Input for the unified design, `_forcing-design.md`
(2026-10-18).

States checked: yelmox dev d60e752, yelmo dev bd9ef661, FesmData (local
checkout: `shared/domains.toml`, dataset READMEs).

## Goal

Two Greenland paleo setups, analogous to `scripts/ant-paleo/`
([Antarctic paleo](experiment-ant-paleo.md)), on the **GreenlandPaleo**
domain (GRL-PAL), with Ellesmere Island and the ice thickness at the domain
borders prescribed in time:

- **A. Glacial cycles**: PD spin-up, then -130 kyr to present with the climate
  blended between snapshots by a glacial index (as ant-paleo).
- **B. Deglaciation after Tabone et al. (2024)**: PD spin-up, LGM steady state,
  then a transient from the LGM to the present forced by the yearly-to-decadal
  reconstructions of Buizert et al. (2018) and Badgeley et al. (2020). Current
  methods and data (chion, FesmData v2); only the forcing method is reproduced.

Both on snapesm (snapclim is being retired, `_snapclim-retirement-plan.md`,
which drops its `recon` mode).

## Domain: GRL-PAL

- Grid (`FesmData/shared/domains.toml`): x [-1360, 960], y [-3450, -570] km,
  polar stereographic as GRL; GRL is its crop x [-720, 960]. Resolutions 0.5
  to 32 km.
- The extra 640 km in the west hold Ellesmere, Devon and Baffin Bay. The west
  border cuts Baffin Island and the Canadian Arctic Archipelago, so
  Innuitian/Laurentide ice crosses it at the LGM: hence prescribed border ice.
- No yelmox config uses GRL-PAL yet.

### FesmData on GreenlandPaleo

| Available | Missing |
|---|---|
| Topo: BedMachine v6 (v5, v4) | Reference climate: MAR v3.14 only on Greenland |
| Atmosphere: ERA5 1981-2010, 1991-2020 | Buizert2018: remapped to GRL-16KM/32KM only; source 96°W-0°, 58-85°N |
| Ocean: ORAS5 | Badgeley2020: remapped to GRL-8/16/32KM only; source (T31) 86.25°W-3.75°W, 53.8-87.2°N |
| Paleoclimate: PMIP3, PMIP4, CLIMBER-3α, CLIMBER-X | TraCE-21ka ocean (not in v2; v1 `OCN-RECON-TRACE21ka`) |
| IceHistory: GLAC-1D (26-0 ka; nn9927, nn9894 120-0 ka), ICE-6G_C (26-0 ka), ICE-5G, PaleoMIST (80-0 ka), Batchelor2019, PaleoGrIS | |
| SolidEarth: Lloyd2024, Bagge2021 (16 km and coarser) | |
| GHF, sediments, basins (Mouginot, Zwally, Schmidt2025, NEGIS) | |

## What recent yelmox changes give

| Change | Use here |
|---|---|
| snapesm index axes, any tsgen method, `const` index (afd4fe3) | A on snapesm, one or two indices (orbital + millennial, Banderas et al., 2018) |
| time-varying `var_defs` groups (varslice `with_time`, refreshed in `snapesm_refresh_loads`, `extrap` = interpolate and clamp) | B: a reconstruction is a time-varying snapshot; `dx = 0` reads one slice (LGM steady state) |
| varslice online remap (`remap`, `grid_src`, `fill`; 1dcb3d3) | reconstructions or global products read on GRL-PAL without regional products |
| `[bias]` static `dtas`/`fpr` (474f1d3) | one route to a MAR-quality reference on GRL-PAL |
| downscaling to the surface grid (3151308) | coarse climate grid, finer ice/surface grid |
| forcing at `time_rel` (7969dcb), `const_kabp` in kyr BP (ed2a6b5) | constant-time LGM spin-up holds 21 ka climate and insolation |
| `init_method = "recon"`, `recon_time` (22a5fa5) | LGM spin-up cold-starts from GLAC-1D/ICE-6G at 21 ka |

## Reference climate on GRL-PAL

MAR covers only the Greenland crop. Options:

1. ERA5 as `pd` reference plus `[bias]` from MAR - ERA5 over Greenland
   (`dtas` 0, `fpr` 1 elsewhere, tapered at the MAR edge). Needs a bias product.
2. A merged MAR/ERA5 product in FesmData (MAR where valid, ERA5 elsewhere,
   blended at the edge).

## Reconstruction forcing in snapesm (B)

Blend: `linear: now = ref + a (hi - lo)`, `ratio: now = ref (1 + a (hi/lo - 1))`.

- **Buizert2018** (absolute monthly `tas`, `pr`, decadal, 22 ka to 2005 CE,
  with a 1981-2010 baseline): `lo` = baseline (static), `hi` = reconstruction
  (time-varying), index `const` = 1. Fits as is.
- **Badgeley2020** (annual `tas` anomaly and `pr` fraction relative to
  1850-2000 CE, 50-yr means, 20 ka to present): already an anomaly; needs a
  zero/one `lo` snapshot or a snapshot option "is an anomaly".
- Combining Buizert temperature with Badgeley precipitation: two snapshots,
  each bound to its own field.
- LGM steady state: the same groups with `dx = 0` (one slice at `x0`), or a
  constant-time run at -21 ka (forcing evaluated at `time_rel`).

To check:

1. Buizert2018 layout (lon, lat, month, time) vs the varslice sub-annual axis
   (`time = 1, x0, x1, dx, 12`): months may have to be folded into the time
   records, in the product or in the reader.
2. I/O per update: two bracketing slices of a 1.5 GB (GRL-16KM) file;
   whether varslice caches them until the bracket changes.
3. Coverage west of 96°W (Buizert) and 86°W (Badgeley): fill (nearest), blend
   to a PMIP/index-scaled anomaly, or taper to the reference.
4. Regional products on GRL-PAL vs online remap from the source grids.
5. Ocean: no TraCE in v2. Options: `fraction` of the atmospheric anomaly,
   CLIMBER-3α/-X ocean snapshots, or a new TraCE product. snapesm needs the
   annual / no-salinity ocean of the retirement plan (step 2) for the
   snapshot oceans.

## Prescribed ice at the domain borders

### In yelmo now

- `bnd%mask_ice`: `MASK_ICE_NONE` (0), `MASK_ICE_FIXED` (1), `MASK_ICE_DYNAMIC`
  (2); fixed cells get `H_ice = bnd%H_ice_ref`
  (`yelmo/src/physics/mass_conservation.f90:529`).
- Set by region expressions in the masks group (`mask_ice_dynamic`,
  `mask_ice_fixed`) and on the border by `yelmo.mask_border`
  (`auto`/`none`/`fixed`/`dynamic`; `auto` has no GRL-PAL case).
- Relaxation zones: `relax`, `relax_tau` -> `bnd%tau_relax` (`topo_rel = -1`),
  towards `H_ice_ref`.

### Decided (2026-10-18)

No separate field: `H_ice_ref` is the prescribed/reference thickness and may
vary in time. Nothing in yelmo assumes it constant; its uses:

| Where | Use | Time-varying value |
|---|---|---|
| yelmo `mass_conservation.f90:529` | `MASK_ICE_FIXED`: `H_ice = H_ice_ref` | the purpose |
| yelmo `yelmo_topography.f90:320/332` | relaxation target (`topo_rel_field = "H_ref"`); `topo_rel` 1/2 masks | fine (PD spin-ups only) |
| yelmox `kryos_startup.f90:210` | isostasy reference, copied once at cold start | fine |
| yelmox `kryos_startup.f90:274-294` | `recon` cold start overwrites it with the reconstruction | already not PD |
| yelmox `kryos_coupling.f90:448` | `smb_simple` target mask, refreshed each step | intended: the mask follows `H_ice_ref` |
| yelmo restart | saved and read | fine |

The optimisation targets `dta%pd%H_ice`, not `H_ice_ref` (the earlier
version of this plan was wrong there).

Proposal:

1. yelmo masks group: `border_width = N M`, a strip of N fixed cells
   (`MASK_ICE_FIXED`) and M relaxation cells inside it (`tau_relax` =
   `relax_tau`, `topo_rel = -1`), measured inward from the edge of the
   `mask_ice_dynamic` region; outside the region `MASK_ICE_NONE`. Applied
   after the expressions (`mask_ice_dynamic`, `mask_ice_fixed`, `relax`),
   which stay for hand-made zones. Covers `mask_border = "fixed"` (`1 0` on
   the domain border).
   - GRL-PAL: `mask_ice_dynamic = "all"`, the strip is the domain border.
   - NH with only the GRL-PAL area active: `mask_ice_dynamic` = an
     expression for that area, from a layer of the regions: a FesmData layer
     (e.g. the footprint of GRL-PAL on the NH grids) or an own mask file
     (`path_<layer>` in the regions group, remapped by nearest neighbour onto
     the grid if needed).
2. yelmox `[ice_bnd]`: a varslice group of an IceHistory product
   (time-interpolated on the forcing clock, `extrap`), PaleoMIST by default,
   any file; writes `H_ice_ref` in the fixed and relaxation cells only (the
   interior keeps its value).
3. Outside the time range of the file (PaleoMIST: 80-0 ka; GLAC-1D nn9927:
   120-0 ka): a scaled thickness, `H = f(index) × H_slice` (e.g. the LGM slice
   times the glacial index), with a `method` key to iterate on the scaling.

Consequences:

- The fixed ice loads FastIsostasy; the bed under it comes from our isostasy,
  not from the reconstruction.
- Sea level: `[barysealevel] method = "file"`; `fastiso` would add the border
  ice to the barystatic sea level.
- Regional output and error metrics exclude the border zone (`mask_rmse`,
  regions).

## Setups

Common: GRL-PAL (16 or 32 km), BedMachine v6, chion, the reference climate
above, ORAS5 ocean reference, FastIsostasy with Lloyd2024 or Bagge2021, sea
level from file, regions separating Greenland, Ellesmere and the border zone.
Scripts in `scripts/grl-paleo/` (or `grl-lgp/` and `grl-deglac/`), par files
`yelmox/yelmox_GreenlandPaleo_*.nml`, docs page `experiment-grl-paleo.md`.

### A. Glacial cycles

1. `run_spinup.sh`: PD, optimise `cb_ref` and `tf_corr` (`L21-points`),
   relaxation early on (as ant-paleo).
2. `run_lgp.sh`: -130 kyr to 0 (or +10), relative time, snapesm `pd`,
   `picontrol`, `lgm` (PMIP3, `clim_sigma` 50 km); index `at`/`ap` from a
   glacial index (e.g. `alpha_combined_125kyr_interp.dat`), optionally `bt`
   from `beta_ngrip_millennial.dat` on a weak/strong LGM axis; transient
   insolation; elsa isochrones (compare with MacGregor 2025); border ice from
   `[ice_bnd]`.

### B. Deglaciation (Tabone et al., 2024)

1. `run_pd_spinup.sh`: as A.1.
2. `run_lgm_spinup.sh`: cold start `init_method = "recon"` (GLAC-1D or ICE-6G
   at 21 ka), constant time at -21 ka: one reconstruction slice, LGM sea level
   (`input/lgm-120m.dat`), border ice at 21 ka; to equilibrium incl.
   isostasy.
3. `run_deglac.sh`: -21/-22 ka to 0, time-varying reconstruction snapshot,
   `input/sl_deglaciation.dat`, transient border ice.

Validation with the metrics of `GrIS_deglaciation_yelmox` (PaleoGrIS extent
and isochrones, Vinther 2009 elevations, Schumacher 2018 GPS uplift).

## Work, in order

1. FesmData: reference climate (or bias) on GRL-PAL; Buizert2018/Badgeley2020
   on GRL-PAL or decision for online remap; western fill; layout fix for
   varslice if needed.
2. yelmo: `border_width` in the masks group.
3. yelmox: `[ice_bnd]`; base `yelmox_GreenlandPaleo.nml` (domain, regions,
   isostasy, masks).
4. yelmox: snapshot axes with `lo = ""` (Badgeley), regional taper modifier;
   see `_forcing-design.md`.
5. A: par files, scripts, docs.
6. B: par files, scripts, docs.

## Decisions (2026-10-18)

1. Border ice: any IceHistory file, PaleoMIST by default; before the file
   starts, a scaled thickness from an index (`[ice_bnd] method`).
2. No `H_ice_fix`: `[ice_bnd]` writes `H_ice_ref` (section "Prescribed ice").
3. Strip of cells, `border_width = N M` (fixed, relaxation), from the edge of
   the dynamic region; own masks via expressions or a regions layer from a
   file; works on GRL-PAL and on NH with only the GRL-PAL area active.
4. Reference climate: ERA5 + MAR - ERA5 bias via `[bias]`. FesmData:
   a bias product on GRL-PAL (`dtas`, `fpr`), tapered at the MAR edge.
5. Reconstructions outside their coverage: regional taper, a climate modifier
   `now = ref + w (now - ref)` (atmosphere), `w` from the `f_valid` of the
   Buizert/Badgeley products (smoothed), no mixing of datasets.
6. Badgeley: snapshot axis with `lo = ""`.
7. GRL-PAL-16KM first, resolution configurable (8 km later).
8. Tabone et al. (2024): temperature and precipitation reconstructions chosen
   separately (one snapshot per field); `cb_ref` optimisation optional.
9. Ocean for B: all three, selectable (`ocn_fraction`, CLIMBER ocean
   snapshots, TraCE-21ka product); see `_forcing-design.md`.
10. One script directory, `scripts/grl-paleo/`.

Still to check (data/code, section "Reconstruction forcing"): the Buizert
month/time layout vs varslice, I/O per update, Buizert/Badgeley on GRL-PAL
(new products or online remap from the GRL grids).
