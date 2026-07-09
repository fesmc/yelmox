---
title: "Program flavors"
---

YelmoX ships several **driver programs** ("flavors"), each a `program` that wires
Yelmo together with a different set of forcing/coupling components. All of the
modern (multigrid) flavors are built on the shared `ice_domain` type and the
`step_*` coupling primitives in [`libs/yelmox_domain.f90`](https://github.com/fesmc/yelmox/blob/main/libs/yelmox_domain.f90);
they differ in **which components are active** and in **how the per-step coupling
sequence is assembled**.

See [Multigrid coupling](multigrid.md) for the design of the shared `ice_domain`
core that these drivers reuse.

## The flavors

| Flavor | Build | Climate / SMB | Ocean | Distinctive feature |
|---|---|---|---|---|
| [`yelmox`](flavor-yelmox.md) | `make yelmox` | snapclim + smbpal | snapclim | Single domain; canonical driver. Transient time-series forcing (`tsgen`). |
| [`yelmox_bipolar`](flavor-bipolar.md) | `make yelmox_bipolar` | snapclim + smbpal (×2) | snapclim + shared OBM | Two hemispheres, shared sea level + Ocean Box Model. |
| [`yelmox_esm`](flavor-esm.md) | `make yelmox_esm` | ESM forcing (`libs/esm.f90`) | ESM | Earth-System-Model forcing replaces snapclim; calendar timeline. |
| [`yelmox_rembo`](flavor-rembo.md) | `make yelmox_rembo` | REMBOv1 | snapclim | REMBO energy/moisture-balance atmosphere + SMB. |

Legacy single-grid originals live under each flavor's `legacy/` folder and build
with `make <flavor>-legacy`. Flavors not yet ported to the multigrid driver
(`yelmox_ismip6`, `yelmox_nahosmip`, `yelmox_rtip`) are legacy-only.

## Shared coupling primitives

Every modern flavor advances the model by calling these primitives (from
`yelmox_domain`), in a flavor-specific order:

- `step_optimize` — spinup relaxation + basal-friction / thermal-forcing tuning.
- `step_isostasy` — bedrock/sea-level (FastIsostasy), against the shared barystatic sea level (`bsl`).
- `step_icesheet` — assemble the Yelmo boundary state (couplers) and run `yelmo_update`.
- `refresh_htopo` — mirror the prognostic geometry into the hi-res reference hub.
- `step_climate` — climate on `grid_clim` + SMB on `grid_smb` (snapclim/smbpal).
- `step_marine_shelf` — sub-shelf melt on `grid_mshlf`.

The single-domain [`yelmox`](flavor-yelmox.md) wraps these in `yelmox_step`; the
other drivers call the primitives directly so they can interleave extra steps
(a second domain, an ocean box model, an ESM/REMBO climate step).

## Initialization ordering (applies to all flavors)

On cold start, the Yelmo applied mass-balance diagnostics (`smb`, `bmb`, `fmb`)
are populated at the initial time so the first output snapshot reflects the
coupled boundary forcing rather than zeros. This is handled inside Yelmo: when
the topography solver runs without advancing the ice (`pc_step="none"` at init,
or any `topo_fixed` step), it diagnoses the applied mass balance from the current
boundary forcing (`calc_ytopo_mb_diagnostic`). `mb_net` remains zero at `t=0`
(nothing is applied), which is expected.
