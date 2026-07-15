# Changelog

All notable changes to YelmoX are recorded here. Each version corresponds to an
annotated git tag. Dates are release (tag) dates.

## [v2.3] - 2026-07-15

### Added
- Multigrid driver framework: all flavors (`yelmox`, `yelmox_bipolar`,
  `yelmox_esm`, `yelmox_rembo`) rebuilt on the shared `ice_domain` type and
  `step_*` coupling primitives in `libs/yelmox_domain.f90`; per-flavor
  documentation pages.
- Transient time-series forcing (`tsgen`/`tsforcing_class`): single forcing value
  mapped onto per-channel anomalies, forcing-increment (`Δf`) restarts, kill
  switch, and 1D diagnostics; wired into the snapclim, ESM, and REMBO drivers.
- `snapesm` climate backend: backend-agnostic `dom%clim` adapter
  (`yelmox_climate`) so an ESM/varslice climate can replace snapclim.
- `yelmox_esm`: ISMIP7 spin-ups (Greenland + Antarctica), TIPMIP Greenland
  stabilisations, and 1pctCO2 forcing-only scaffolds; annual-mean transient ocean
  via `varslice_nsub`.
- Namelist database in the compact fesm-utils:dev `varslice` format, with a
  database-namelists documentation page.

### Changed
- `yelmox_esm`: `esm` module/source renamed to `esm_forcing` (avoids ifx name
  clash); forcing loaded through shared `varslice_init_nml` with
  `{gcm}`/`{experiment}` substitutions.
- Sync with yelmo:dev: `ytrc` tracer subsystem, `var_io` tables, and default
  parameters.
- `yelmox_rembo`: driver-contained routines take explicit arguments instead of
  relying on host association, matching `yelmox_esm`/`yelmox_bipolar`.
- Output layout: drop grid suffix, add `yelmo_sm.nc`, htopo mask-load flags.

### Removed
- Retired the `yelmox_ismip6`, `yelmox_nahosmip`, and `yelmox_rtip` flavors to
  `retired/`. They still compile (`make yelmox_ismip6` / `yelmox_nahosmip` /
  `yelmox_rtip`, each printing a retirement notice); prefer `yelmox` or
  `yelmox_esm` for new work.
- Dropped pre-configme / pre-runme tooling.

### Fixed
- `var_io` tables (`input/yelmo-variables-{ydata,ytrc}.md`) synced to the
  yelmo:dev isochrone dimension rename (`age_iso`→`time_iso`,
  `pd_age_iso`→`pd_time_iso`). The stale copies crashed restart writing
  (`nf90_inq_dimid`), which affected every run.

## [v2.2.2] - 2026-06-24
- yelmox: added bsl ts writing.

## [v2.2.1] - 2026-06-24
- Config and small bug fixes.

## [v2.2] - 2026-06-18
- Yelmo default parameters added. `yhyd` section (with bug fix to bucket units).
  Added yelmo-config tool for inspecting and comparing parameters. yelmox: new
  capabilities to support ISMIP7 activities.

## [v2.1.3] - 2026-06-15
- Added support for ISMIP7 simulations GrIS+AIS.

## [v2.1.2] - 2026-06-15
- Added `with_isostasy` parameter for yelmox - now complete.

## [v2.1.1] - 2026-06-15
- Added `with_isostasy` parameter for yelmox.

## [v2.1] - 2026-06-11
- Conversion to FastHydrology for yelmo basal hydrology.

## [v2.0.6] - 2026-06-02
- Added missing yelmo variables in io tables.

## [v2.0.5] - 2026-06-01
- Added missing isos parameter, removed obsolete.

## [v2.0.4] - 2026-05-30
- yelmox/yelmo work with configme v0.6.6+ and runme v0.5.8+.

## [v2.0.3] - 2026-05-30
- yelmox/yelmo work with configme v0.6.2+ and runme v0.5.8+.

## [v2.0.2] - 2026-05-24
- Tagged version consistent with yelmo:v2.0.2, with central runme and configme setup.
