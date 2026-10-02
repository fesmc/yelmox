---
title: "Program flavors"
---

YelmoX ships several **driver programs** ("flavors"), each a `program` that wires
Yelmo together with a different set of forcing/coupling components. All of the
modern (multigrid) flavors are built on the shared `kryos_domain` type and the
`step_*` coupling primitives of the Kryos modules (`libs/kryos*.f90`);
they differ in **which components are active** and in **how the per-step coupling
sequence is assembled**.

See [Multigrid coupling](multigrid.md) for the design of the shared `kryos_domain`
core that these drivers reuse.

## The flavors

| Flavor | Build | Climate / SMB | Ocean | Distinctive feature |
|---|---|---|---|---|
| [`yelmox`](flavor-yelmox.md) | `make yelmox` | `[coupling] climate`: snapclim, snapesm or [esm](flavor-esm.md); smbpal, smb_simple or the climate's smb | the climate (profiles or shelf base) | Single domain; canonical driver. Transient time-series forcing (`tsgen`). |
| [`yelmox_bipolar`](flavor-bipolar.md) | `make yelmox_bipolar` | snapclim + smbpal (×2) | snapclim + shared OBM | Two hemispheres, shared sea level + Ocean Box Model. |
| [`yelmox_rembo`](flavor-rembo.md) | `make yelmox_rembo` | REMBOv1 | snapclim | REMBO energy/moisture-balance atmosphere + SMB. |

## Shared coupling primitives

Every modern flavor advances the model by calling these primitives (from
`kryos_coupling`), in a flavor-specific order:

- `step_spinup_tuning` — spinup relaxation + basal-friction / thermal-forcing tuning.
- `step_isostasy` — bedrock/sea-level (FastIsostasy), against the shared barystatic sea level (`bsl`).
- `couple_to_yelmo` — assemble the Yelmo boundary state from the component outputs (incl. the climate's subglacial discharge, when supplied).
- `step_icesheet` — run `yelmo_update`.
- `refresh_hub` — the hub's current geometry from the models (a mirror of Yelmo on its grid; hi-res reference + Yelmo's anomalies on a finer hub).
- `step_climate` — climate on `grid_clim` from the backend chosen by `[coupling] climate` (`snapclim`, `snapesm` or `esm`), with the transient forcing, on the `dt_clim` cadence.
- `step_smb` — surface mass balance on `grid_smb` (`smb_method`: smbpal, smb_simple, or the climate's own, `climate`).
- `step_marine_shelf` — sub-shelf melt on `grid_mshlf`, from the climate's ocean as depth profiles or at the shelf base.

Every driver writes the sequence out in its time loop, so the coupling order can
be read directly from the program; drivers with extra steps (a second domain, an
ocean box model, an ESM/REMBO climate step) interleave them there.

## Initialization ordering (applies to all flavors)

On cold start, the Yelmo applied mass-balance diagnostics (`smb`, `bmb`, `fmb`)
are populated at the initial time so the first output snapshot reflects the
coupled boundary forcing rather than zeros. This is handled inside Yelmo: when
the topography solver runs without advancing the ice (`pc_step="none"` at init,
or any `topo_fixed` step), it diagnoses the applied mass balance from the current
boundary forcing (`calc_ytopo_mb_diagnostic`). `mb_net` remains zero at `t=0`
(nothing is applied), which is expected.
