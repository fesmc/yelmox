---
title: "yelmox_bipolar"
aliases:
  - flavor-bipolar.html
---

The two-domain variant of [`yelmox`](yelmox.md): it runs a northern and a southern
ice-sheet domain together, coupled through a **shared barystatic sea level** and a **shared Ocean
Box Model (OBM)** that exchanges freshwater flux and ocean temperature between the
hemispheres. Each domain is a full `kryos_domain`, configured
as in `yelmox` with the hemisphere suffix on its groups (see
[Configuration](yelmox.md#configuration)); the driver interleaves the `step_*` primitives
across both domains plus the OBM.

- **Program:** `yelmox_bipolar/yelmox_bipolar.f90` + `yelmox_bipolar/obm_coupling.f90` + `libs/kryos*.f90`.
- **Build:** `make yelmox_bipolar` (links the OBM stack, `$(obm_libs)`).
- **Config:** `yelmox_bipolar/yelmox_bipolar_Bipolar.nml`: north = Greenland
  (hub and ice `GRL-16KM`, the other components `GRL-32KM`), south = Antarctica
  (hub and shelf `ANT-16KM`, the other components `ANT-32KM`, with `relax` and
  `opt`); `[ctrl]` 0 to 40 yr, `dtt = 10` yr, the box model off (`active_obm`,
  `atm2obm`, `ism2obm` = False). `[ctrl]`, `[barysealevel]` and the Yelmo
  physics groups (`ydyn`, `ytopo`, ...) are shared (no suffix).
- **Tests:** `yelmox_bipolar/run_tests.sh` (local runs under `tmp/`): coupling
  off; the box model on (`atm2obm` + the `obm_north`/`obm_south` series and
  terms in place of the `ocean` term);
  with `RUN_FWF=1` also `ism2obm` (needs the hydrographic masks).

## What's distinct

- **Two `kryos_domain`s** (`dom_north`, `dom_south`), each set up via `setup_domain`
  → `domain_startup`. Either can be individually deactivated (`active_north` /
  `active_south`); with neither, the ocean box model runs alone (`active_obm`,
  forced by `atm2obm` or not; `ism2obm` off, no terms on the box series).
- **Shared `bsl`** — one barystatic sea level for the run, restored once at startup
  (`bsl_startup`) and written to the run restart bundle.
- **Shared OBM** (`obm`) — an ocean box model stepped once per timestep, configured
  via `obm_ctl_load` and coupled to the domains through `obm_coupling.f90`
  (`coupling_obm2ism`, `obm_coupling_init`, `obm_exchange`). The OBM writes its own
  1D output, `<obm_name>.nc` (`[ctrl] obm_name`, `nautilus` or `stommel`), on the
  `[tm_1D]` cadence; its restart (`obm_restart.nc`) goes into the run restart
  bundle next to `bsl_restart.nc`.
- **One restart bundle per run** — `restart-<kyr>-kyr/` holds the shared
  `bsl_restart.nc`, `obm_restart.nc` and the series (`series_<name>_restart.nc`),
  and each domain in a subfolder named after it (`Greenland/`, `Antarctica/`).
  `[ctrl] restart` is the bundle to start from. The output of each domain goes
  to the same subfolder of the run folder.

Climate/SMB per domain is **snapshots + chion** (`[chion_north] model = "itm"`,
`[chion_south] model = "pdd"`), as in the single-domain driver;
the OBM's contribution is folded into the ocean forcing — `obm_exchange` sets
the derived series `obm_north`/`obm_south` (the box temperatures) and re-applies
each domain's ocean terms (e.g. the `replace` terms `term_obm_north`/`_south`:
the box temperature into `to_ann`) before the marine-shelf step reads it.

## Stepping order

Main loop (per timestep):

```fortran
call bsl_update(bsl, ts%time_rel)              ! shared sea level, once
call indices_update(idx, ts%time_rel, var=fvar) ! shared series, once (fvar: ice volume of the active domains [Gt])

if (active_north) then                         ! relaxation + optimization + isostasy, per domain
    call step_relax(dom_north, ts)
    call step_optimize(dom_north, ts)
    call step_isostasy(dom_north, ts, bsl)
end if
(same for dom_south)

if (oc%active_obm) call obm_update(obox, real(dtt,wp), oc%obm_name)   ! ocean box model, one step

if (active_north) then                         ! ice sheet, hub, climate + smb, per domain
    call couple_to_yelmo(dom_north)
    call step_ice_bnd(dom_north, ts, idx)
    call step_icesheet(dom_north, ts)
    call couple_yelmo_to_htopo(dom_north)
    call step_climate(dom_north, ts, idx)
    call step_surface(dom_north, ts)
end if
(same for dom_south)

call obm_exchange(oc, obox, idx, dom_north, dom_south, ...)  ! atm->obm, ism->obm freshwater,
                                                            ! hysteresis forcing, obm->ism (the obm_* series, ocean terms)

if (active_north) call step_shelf(dom_north, ts)  ! reads the obm-updated to_ann
if (active_south) call step_shelf(dom_south, ts)
```

Key ordering points:

- **Isostasy for both domains runs before the OBM step**, so the OBM sees a
  consistent geometry.
- **`obm_update` uses the previous step's** atmospheric/freshwater forcing (a
  one-step lag), then `obm_exchange` distributes the fresh OBM state back to the
  domains before the marine-shelf melt is computed.
- The per-domain primitives are the same as in the single-domain `yelmox`; only
  the relaxation + optimization + isostasy part is split off so it runs before the OBM step.

## Forcing

One global series drives the box model and, with regional gains, the domains
(`yelmox_bipolar_Bipolar.nml`):

- **Series** (`[forcing] series`, `[series_<name>]`) are shared by both
  domains and the OBM: advanced once per step, written to each domain's 1D
  file, restarted from the run bundle (`series_<name>_restart.nc`), with the
  kill switch and `restart_every_df` bundles as in `yelmox`. Here `alpha`
  (the orbital record r2020, normalised, of the north's term `lgm`), `gmt` (global mean air
  temperature anomaly [K], `const` 0; a file series instead for a transient)
  and `dTo` (derived, the north's ocean anomaly). The reference and the
  terms of each domain are `[forcing_north]` and `[forcing_south]`.
- **`atm2obm`** (`[ctrl]`): the box air temperatures from `[obm_forcing]`,
  `theta_i = theta_i_init + ampl_i × gain × index` for north, tropics and south
  (`ampl_i`: `[nautilus] thermal_ampl_*`, Zickfeld et al., 2004) and the vapour
  fluxes `phin`, `phit` linear in it. It needs no ice-sheet domain; the series
  must be in the registry (checked at init).
- **Domains**: the uniform terms `regional_north`/`regional_south` on `gmt`
  with `f_ta` the regional amplification (north 1.07, south 0.86): explicit in
  the par file, zero while `gmt` = 0. North: `pd_grl` plus the term `lgm`
  (picontrol -> lgm) and the term `ocean` (the derived series `dTo`); south:
  `pd_ant` only.
- **`obm2ism`**: per domain, a series derived from the box
  (`[series_obm_north]`, `[series_obm_south]` of the library: `method =
  "derived"`, `from = "obm_north"`/`"obm_south"`, the box temperature, or its
  anomaly to the initial value with `anomaly = True`) with a term on `to`
  (`[term_obm_north]`/`[term_obm_south]`: `replace`, `to_ann` = the box
  temperature; a `uniform` term with the anomaly series for the reference +
  the box anomaly). `obm_exchange` sets the series and re-applies the ocean
  terms after both climate steps. At the start the box model is initialized
  (cold start or restart) and the series set before the domains, so their
  startup climate, shelf and ice state already see the box. It needs `active_obm`; the series replace `dTo`
  (`forcing_north.terms==lgm,regional_north,obm_north`, see
  `run_tests.sh`).
- **`ism2obm`**: freshwater flux of each domain into the boxes (`fn`, `fs`),
  per hemisphere with `[ctrl] couple_fwf_north`/`couple_fwf_south`, from the
  Yelmo mass balance (`fwf_definition`: `dVdt`, `mb`, `ablation`) within the
  hydrographic mask (`hydro_mask_north`/`hydro_mask_south`, on the Yelmo grid).
- **Hysteresis** (`nautilus` only): `[nautilus] hyster_on` nudges one box-model
  control (`hyster_forcing`: `phit`, `phin`, `fn`, `fs`, `global_melt`) along
  `hyster_forcing_method` (`triangular`, `sin`, `noise`, `linear`) in
  `obm_exchange`.

`obm_restart.nc` (run bundle) holds the box state, the fluxes and the box air
temperatures; the `*_init` values stay those of the par file, the reference of
`atm2obm` and of the box anomaly (`anomaly = True`).
