---
title: "Programs and climate backends"
---

YelmoX has two **programs**, both built on the shared `kryos_domain` type and the
coupling primitives of the Kryos modules (`libs/kryos*.f90`):

| Program | Build | Domains | Climate | Ocean |
|---|---|---|---|---|
| [`yelmox`](flavor-yelmox.md) | `make yelmox` (`rembo=1` for REMBO) | one | any backend (`[coupling] climate`) | the climate (depth profiles or shelf base) |
| [`yelmox_bipolar`](flavor-bipolar.md) | `make yelmox_bipolar` | north + south | snapclim (×2) | snapclim + shared Ocean Box Model |

The forcing of a `yelmox` run is set at runtime by its **climate backend**,
`[coupling] climate`:

| Backend | Supplies | Page |
|---|---|---|
| `snapclim` | atmosphere + ocean profiles, from climate snapshots blended by indices | [snapclim and snapesm](climate-snap.md) |
| `snapesm` | as snapclim, configured as one blend model over any number of snapshots | [snapclim and snapesm](climate-snap.md) |
| `esm` | reference climatology + Earth-system-model anomalies; atmosphere or smb, ocean at the shelf base, subglacial discharge | [ESM forcing](flavor-esm.md) |
| `rembo` | REMBOv1 atmosphere + smb; ocean from snapclim | [REMBO climate](flavor-rembo.md) |

See [Multigrid coupling](multigrid.md) for the architecture shared by both programs.

## Shared coupling primitives

Both programs advance a domain with these primitives (from `kryos_coupling`):

- `step_spinup_tuning` — spinup relaxation + basal-friction / thermal-forcing tuning.
- `step_isostasy` — bedrock/sea-level (FastIsostasy), against the shared barystatic sea level (`bsl`).
- `couple_to_yelmo` — assemble the Yelmo boundary state from the component outputs (incl. the climate's subglacial discharge, when supplied).
- `step_icesheet` — run `yelmo_update`.
- `couple_yelmo_to_htopo` — the hub's current geometry from the models (a mirror of Yelmo on its grid; hi-res reference + Yelmo's anomalies on a finer hub).
- `step_climate` — climate on `grid_clim` from the backend, with the transient forcing, on the `dt_clim` cadence.
- `step_smb` — surface mass balance on `grid_smb` (`smb_method`: smbpal, smb_simple, or the climate's own, `climate`).
- `step_marine_shelf` — sub-shelf melt on `grid_mshlf`, from the climate's ocean as depth profiles or at the shelf base.

Every driver writes the sequence out in its time loop, so the coupling order can
be read directly from the program; `yelmox_bipolar` interleaves the second domain
and the ocean box model there.

## Initialization ordering

On cold start, the Yelmo applied mass-balance diagnostics (`smb`, `bmb`, `fmb`)
are populated at the initial time so the first output snapshot reflects the
coupled boundary forcing rather than zeros. This is handled inside Yelmo: when
the topography solver runs without advancing the ice (`pc_step="none"` at init,
or any `topo_fixed` step), it diagnoses the applied mass balance from the current
boundary forcing (`calc_ytopo_mb_diagnostic`). `mb_net` remains zero at `t=0`
(nothing is applied), which is expected.
