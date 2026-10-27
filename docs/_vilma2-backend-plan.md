# Global sea level and regional isostasy (VILMA2): plan

Status: design, revised 2026-10-27 after an independent review (first
version 2026-10-23). The leading `_` keeps quarto from rendering it.

States checked: yelmox dev f9fc3d0, FastIsostasy dev 22a458f (yelmox
submodule), vilma2 main 60a85df, yelmo dev, climber-x 01cfd63, FesmData
(local).

## Goal

Add VILMA2 (global, spherical, gravitationally self-consistent sea level)
to yelmox next to FastIsostasy (regional bedrock and regional geoid), so
that each can run alone or both together.

FastIsostasy and VILMA2 are not two backends for one component. VILMA2 is
a global Earth and sea-level model, FastIsostasy a regional Earth model.
yelmox therefore gets two layers:

- **Global sea level**, one per run, always present. It holds the global
  ice distribution and determines the sea level from it.
- **Regional isostasy**, one per domain, optional. It computes the bedrock
  inside the domain (and, without VILMA2, the regional sea-surface shape),
  given the global sea level.

```
[sealevel]
    method = "vilma2"         ! "const" | "file" | "barystatic" | "vilma2"

[isostasy]                    ! per domain ([isostasy<sfx>])
    method = "fastisostasy"   ! "none" | "fastisostasy" | "global"
```

The defaults are the most complete physics: VILMA2 for the global sea level
and far field, FastIsostasy in the domains. Today's yelmox (one domain,
FastIsostasy, `fastiso` sea level) is the special case `barystatic` +
`fastisostasy`: reproduced bit for bit by phase 1a, closely but not bit for
bit from phase 1b on (see Phases). Every par file sets `[sealevel] method`
explicitly, so a default never changes a configured experiment.

VILMA2 is a default build dependency (`vilma=1`); `vilma=0` builds a stub
for machines without it.

## Decisions

2026-10-23 (still valid):

- Grids are independent; the coupler (fesm-utils `coords`) remaps.
- VILMA2 is extended to return bed uplift and sea-surface change
  separately, and the barystatic sea level.
- Ice outside the domains: a fixed reference (`ice_ext = "ref"`) or a
  global ice history (`"file"`).
- Global reference topography and the 3D viscosity file are namelist
  options with FesmData defaults.
- `dzbdt_corr` belongs to FastIsostasy; it is ignored, with a warning when
  non-zero, in domains that do not run FastIsostasy.
- The full `&vilma` group is configurable from the par file. Its defaults
  are `input/vilma_defaults.nml`, a synchronized copy of VILMA2's (the
  chion/elsa convention).

2026-10-27:

- Two layers instead of two backends: `kryos_sealevel` (run) and
  `kryos_isostasy` (domain) are the only modules that touch VILMA2,
  FastIsostasy or `bsl_class`. VILMA2 and FastIsostasy can run together.
- Sea level comes only from the global ice distribution and the Earth's
  response to it. Where a cell's ice came from (Yelmo or external data)
  does not enter.
- A domain's footprint (the cells whose ice is Yelmo's) is where Yelmo's
  ice is dynamic (`mask_ice == MASK_ICE_DYNAMIC`); every other cell, also
  inside the Yelmo domain, takes the external ice. (See open question 1 on
  `MASK_ICE_FIXED`.)
- Grounded/floating and the coastline inside a Yelmo domain are decided at
  Yelmo's resolution (`f_grnd`), not by the global model on its coarse
  grid. VILMA2 is extended to take them from the host.
- FastIsostasy + VILMA2 includes the far-field bed response, by linear
  superposition in VILMA2.
- `fastiso` and `mixed` are replaced by `barystatic`. Ice outside the
  domains is always a field (reference or history), never a sea-level
  series.
- Global grid `[sealevel] grid`, default `GLOBAL-1DEG` (FesmData lon-lat),
  for every method. VILMA2 maps it to its own Gauss grid internally.
- Defaults are the most complete physics: `[sealevel] method = "vilma2"`,
  `[isostasy] method = "fastisostasy"`, `ice_ext = "file"`.
- VILMA2 runs as in CLIMBER-X (`nml/geo_par.nml`, `&vilma`): `lmax = 170`,
  PREM incompressible (80 km elastic lithosphere, 9 Maxwell mantle layers,
  "Clemens" 1-D reference viscosity), 3D viscosity on, Bagge et al. (2021)
  structure 6 (`v_1.0_s16`). CLIMBER-X's `input/vilma2/visc3d_Bagge2021.nc`
  is identical to it (cdo: maximum difference 0).
- FastIsostasy keeps its own flotation (it does not take Yelmo's `f_grnd`).

## Current coupling (FastIsostasy)

Per domain, `dom%isos` is a `fastisostasy::isos_class` on `grid_isos`. The
driver owns one `bsl_class`, shared by all domains.

| Call site | What |
|---|---|
| `kryos.f90:348` | `isos_init(dom%isos, ..., nx_i, ny_i, dx_isos, dy_isos, ...)` |
| `kryos_startup.f90:146,218` `domain_init_isostasy` | `isos_init_ref` (Yelmo `z_bed_ref`, `H_ice_ref`) + `isos_init_state` (Yelmo `z_bed`, `H_ice`); before `step_ice_bnd` (153), `yelmo_init_state` (179), `domain_init_ice` (182) |
| `kryos_startup.f90:~460` `domain_restart_read` | `isos_init_state` with `use_restart` |
| `kryos_startup.f90:421` | `isos_restart_write` |
| `kryos_startup.f90:42,106` | `bsl_startup`, `run_restart_write` (`bsl_restart.nc`) |
| `kryos_coupling.f90:119` `step_isostasy` | on `dt_isos`: remap `H_ice`, `dzbdt_corr` to `grid_isos`; `isos_update(..., bsl, dwdt_corr)` |
| `kryos_coupling.f90:186` `couple_isostasy_to_yelmo` | `z_bed = z_bed_ref + remap(w + we)`; `z_sl = bsl + remap(z_ss - bsl)` |
| `kryos_coupling.f90:229` `check_isostasy_reference` | `out%z_bed - (w + we)` vs Yelmo `z_bed_ref` on `grid_isos` |
| `kryos_output.f90:260,383` | `isos.nc`, `isos_ts.nc` |
| `yelmox.f90:92,183,197`, `yelmox_bipolar.f90:125,212,225` | `bsl_init`, `bsl_update`, `step_isostasy` |

`[bsl] method` in the par files: `fastiso` everywhere except SRG (`const`)
and GreenlandPaleo (`file`). `mixed` is unused.

Only anomalies cross into Yelmo (`dz_bed`, `dz_ss`), so Yelmo keeps its own
high-resolution reference bed. That contract carries over unchanged.

Facts the design depends on:

- **Bipolar order.** Under `fastiso`, `isos_update` writes the domain's
  contribution into the shared bsl, so north sees south's contribution one
  step late. Under `mixed`, `bsl_update` resets bsl from the file every
  step, so north never sees south's contribution
  (`fastisostasy.f90:843`). The new step order removes both.
- **Masked reference ice.** Yelmo sets `H_ice_ref = dta%pd%H_ice` without a
  mask (`yelmo_ice.f90:1267`), while `H_ice = 0` in `MASK_ICE_NONE` cells.
  FastIsostasy's reference uses the former and its load the latter, so in
  such cells it sees a permanent unloading (and, under `fastiso`, a bsl
  contribution at init). The splice removes this.
- **`H_ice_ref` is not a fixed reference.** `ice_bnd_update` overwrites it
  in fixed and relaxation cells every step (`ice_bnd.f90:133`), and
  `domain_init_recon` sets it to the reconstruction (`init_method =
  "recon"`, LIS and North). FastIsostasy is unaffected only because it
  captures its reference before both and keeps it in `isos_restart.nc`.
- **Init cost of the bsl.** `isos_init_state` calls `isos_update` with
  `advance_bsl` true on a cold start and `V_old = 0`
  (`fastisostasy.f90:544-549, 621-631`), so bsl changes at init whenever
  the initial state differs from the reference.

## Global sea-level layer

`libs/kryos_sealevel.f90`, type `sealevel_class`, owned by the driver. It
replaces `bsl_class` in the drivers. One per run, shared by all domains.

### Densities

The budget uses one pair of densities, the host's (Yelmo's physical
constants: `rho_ice = 910`, `rho_sw = 1028`; FastIsostasy takes its
densities from the same constants). VILMA2 hard-codes 931 and 1000 (the GIA
benchmark convention, `vilma_constants.f90:27`) for load, melt source and
flotation; they become parameters that yelmox sets. bsl is in metres of
seawater at the host's `rho_sw`.

### The global ice field

- **Host region** of domain d: the whole Yelmo domain. Its weight `w_d` on
  the global grid is the conservative remap of ones. Host regions must not
  overlap (checked at init); this also keeps FastIsostasy regions apart.
- **Footprint** `m_d` (on the Yelmo grid): the cells whose ice is Yelmo's,
  `mask_ice == MASK_ICE_DYNAMIC` at the end of init. It is fixed for the
  run; a later change of `mask_ice` there stops the run.
- **Inside the host region**, the ice is assembled on the Yelmo grid:
  `H = m·H_yelmo + (1 - m)·H_ext`, with the external ice interpolated to
  the Yelmo grid. The host-resolved terms (below) are computed there for
  every cell and remapped conservatively. In global cells only partly in
  the host region, the remainder `1 - w_d` holds the external field; both
  parts are then external data, so nothing is counted twice and no coarse
  external anomaly overlaps Yelmo's ice.
- **Outside host regions**: the external ice on the global grid, with the
  global model's own flotation and coastline rule.
- **External ice**:
  - `ice_ext = "ref"`: `H_ref,ext`, fixed;
  - `ice_ext = "file"`: `H_ref,ext + [H_hist(t) - H_hist(t_ref)]`, clipped
    at 0 (the clipped mass is reported). The history enters as an anomaly
    against its own present day, so a mismatch between the history and the
    reference product cannot create a load. After the history's last time
    its last field is held (projections run); before its first time the run
    stops.
- **Reference state**, fixed for the run, with the same assembly: Yelmo's
  present-day data (`dta%pd`: `z_bed`, `H_ice`) in the footprint, the
  external reference elsewhere. `dta%pd` holds exactly what FastIsostasy
  takes as its reference today, but `bnd%H_ice_ref` does not stay equal to
  it (see above). `z_sl_ref = 0` (present-day sea level). `f_grnd_ref`
  comes from Yelmo's own routines on the reference fields (`calc_H_grnd`
  without `f_ice`, then `calc_f_grnd_subgrid_area` or `_linear` as the
  run's `gl_sep`). The reference is rebuilt identically on a restart.

**Host-resolved terms**, computed on the Yelmo grid for every cell of the
host region. `f_g` is Yelmo's `f_grnd` in the footprint, and Yelmo's
flotation rule applied to `H_ext`, `z_bed`, `z_sl` elsewhere:

| Term | Field on the Yelmo grid |
|---|---|
| grounded-ice column `I_g` | `H · f_g` |
| ocean fraction `C` (open ocean and floating ice) | `1 - f_g` (ice-free land has `f_g = 1`) |
| flooded-area fill `F = ζ⁽⁰⁾(C - C⁽⁰⁾)` | `(z_bed_ref - z_sl_ref)·(C - C_ref)` |

`I_g` and `F` are mass fields and are remapped conservatively. `C` and
`w_d` are fractions: conservative too, but kept in [0, 1] and never
rescaled by a global mass factor. `C` is fractional in cells with a coast,
a grounding line or a host-region edge; the budget below only uses
pointwise products and integrals, so that is fine. Products with the
smooth sea-level field (`C·S`) are formed on the global model's grid. The
global integrals use true (spherical) areas; Yelmo's own integrals use
projected `dx·dy`, which differs by up to about 5% on the polar
stereographic grids.

FastIsostasy receives the same assembled `H` on `grid_isos` (conservative
when coarser, bilinear when finer, which is then not mass-exact) and
keeps its own flotation.

**Limitation.** Closed depressions below sea level count as ocean, in
Yelmo's `f_grnd`, VILMA2's `ocean_function` and FastIsostasy alike, so a
deglaciated interior basin adds ocean volume through `F`. A connectivity
test (as CLIMBER-X's `l_connect_ocn`) is a later option.

### One budget for barystatic and vilma2

With `S` the relative-sea-level change (sea surface minus bed), `N` the
sea-surface (geoid) change, `u` the bed change and `Δφ` the uniform shift,
VILMA2 solves `S = N - u + Δφ` with ocean-mass conservation
(`vilma_sle.f90:291-336`, Martinec 2018 eqs 17-20):

    ∫ C·S dA = -(ρ_i/ρ_w) ∫ ΔI_g dA + ∫ F dA

For every method the layer reports:

- `bsl = ∫ C·S dA / ∫ C dA`: the mean relative-sea-level change over the
  ocean (VILMA2's "barystatic rise", the ocean mean of `C·rsl`). Water
  that floods newly submerged area (`F`) does not raise it. A `file`
  series is taken to be bsl in this sense.
- `bsl_ive = -(ρ_i/ρ_w) ∫ ΔI_g dA / ∫ C dA`: the ice-volume equivalent
  (VILMA2's current `se%bsl`), as a diagnostic.
- `dz_ss = N + Δφ`: the sea-surface change.
- `dz_bed = u`: the bed change.

No term depends on where a cell's ice came from: bsl, `Δφ` and `dz_ss` are
set by the global ice field and the Earth's response to it.

**`barystatic`** is this budget with a rigid Earth, except where a regional
isostasy runs. There `u = u_R` and `N = N_R` (FastIsostasy's bed and
regional geoid from its last update, zero outside its regions), so

    Δφ = bsl + ∫ C·(u_R - N_R) dA / ∫ C dA

and each domain gets `dz_ss = N_R + Δφ`, consistent with FastIsostasy's
`z_ss = shift + ref%z_ss + dz_ss` (`sealevel.f90:32`, `ref%z_ss = 0`). The
coastline outside host regions is iterated as in VILMA2, starting from the
previous update's coastline. The state is small: `Δφ`, `bsl` and `C` of
the last update, in `sealevel_restart.nc`, so a restart redo uses the
restored values rather than recomputing them from the post-update state.

**`vilma2`** is VILMA2's sea-level equation with the host-resolved terms
inside host regions.

### Methods

| Method | bsl | Sea surface for domains | Bed for domains |
|---|---|---|---|
| `const` | `bsl_init` | uniform | none |
| `file` | series (`sl_path`) | uniform | none |
| `barystatic` | rigid-Earth budget | uniform shift with regional terms | none |
| `vilma2` | VILMA2 budget | VILMA2 field | VILMA2 field (`global`) or far field (`fastisostasy`) |

Under `const` and `file`, `Δφ = bsl`.

The global ice field is built for every method and written to
`sealevel.nc`. Under `const` and `file` the barystatic budget is also
written, as a diagnostic: it shows how far a prescribed sea level is from
the modelled ice.

## Regional isostasy

`libs/kryos_isostasy.f90`, per domain (`dom%isos`, type
`isos_domain_class`: method, the FastIsostasy state `fi`, and `dz_bed`,
`dz_ss` on the Yelmo grid).

| `[isostasy] method` | `dz_bed` | `dz_ss` |
|---|---|---|
| `none` | 0 | from the global layer |
| `global` (needs `vilma2`) | VILMA2 `u` | VILMA2 |
| `fastisostasy`, global not `vilma2` | `w + we` | `N_R + Δφ` (own regional geoid) |
| `fastisostasy` with `vilma2` | `w + we + u_far` | VILMA2 |

FastIsostasy gets:

- load: the assembled ice field on `grid_isos`;
- sea level: the uniform shift `Δφ` (today's `bsl` argument, read only);
- under `vilma2`: VILMA2's sea-surface anomaly in place of its own
  regional geoid, and the far-field bed `u_far` as an external
  displacement. `u_far` enters FastIsostasy's bed for flotation and the
  ocean column, but not its own solve.

`couple_isostasy_to_yelmo` keeps its formulas,
`z_bed = z_bed_ref + remap(dz_bed)` and `z_sl = Δφ + remap(dz_ss - Δφ)`.
The remaps move into `isostasy_step`, which fills `dom%isos%dz_bed/dz_ss`
on the Yelmo grid (FastIsostasy → Yelmo: `remap_method_smooth` as now;
global → Yelmo: `bilin`).

`check_isostasy_reference` subtracts `u_far` as well, and also compares
FastIsostasy's reference ice with the assembled reference.

### Step order

```fortran
! per domain: relaxation, optimization (as now)
call sealevel_step(sl, doms, ts)       ! once per run: assembly, budget or VILMA2 solve
! per domain:
call isostasy_step(dom, sl, ts)        ! FastIsostasy, or land the global fields
call couple_to_yelmo(dom)
```

The sea level is settled from all domains before any regional solve.
Cadences: `[sealevel] dt` for the global layer (from phase 1b; in phase 1a
the provisional barystatic follows each domain's `dt_isos`); FastIsostasy
keeps its per-domain `dt_isos`. Global fields are held between global
steps.

### Init and restart

Cold start:

1. Driver: `sealevel_init` reads `[sealevel]`, the global grid and the
   external data, and sets the initial shift `Δφ0` (`const`/`file`:
   `bsl(t0)`; `barystatic`/`vilma2`: 0, the reference).
2. Per domain, `domain_startup` as today. Where `domain_init_isostasy` runs
   now: the reference snapshot (`dta%pd`, `f_grnd_ref`) and FastIsostasy's
   reference and initial state, with `Δφ0`. Then `step_ice_bnd`,
   `yelmo_init_state` and `domain_init_ice` as today, on that bed.
3. Driver, after all domains: footprints from the final `mask_ice`, the
   global reference and the `t0` ice. `barystatic`: `Δφ(t0)`. `vilma2`:
   `solid_earth_init` with the reference and host terms, then the spin-up
   (with the channels) and the first fields.
4. Per domain: land the global fields (`dz_ss`, `u_far`, `Δφ`) and
   reconcile Yelmo's topography diagnostics (`calc_ytopo_diagnostic`), as
   `domain_init_isostasy` does now.

VILMA2's spin-up equilibrates the far field and the channels. FastIsostasy
cold-starts with `w = 0`, so the near field starts at rest, as today.

Restart:

1. Driver: read `sealevel_restart.nc` or `vilma_restart.nc`: `Δφ`, bsl,
   `C`, the host terms of the last update (VILMA2's ramp start), VILMA2's
   state and the channel memories.
2. Per domain: `domain_restart_read`; `isos_init_state` with `use_restart`
   gets the restored `Δφ`, sea surface and `u_far`.
3. The first step (`time == time_init`) does not recompute the global
   layer, matching FastIsostasy's restart gate.

In `yelmox_bipolar`, steps 1 and 3 run once in the driver, around the two
domain startups.

## FastIsostasy with VILMA2: the far field

The bed in a domain is the response to the whole global load. With both
models, FastIsostasy takes the near field (its region `R_d`, at its
resolution and with its Earth structure) and VILMA2 the rest.

**Matching load definitions.** FastIsostasy's load is `ρi·ΔHaf` plus
`ρw·(z_ss - z_ss_ref)` on its ocean mask (`sealevel.f90:62-77`). Its ice
term equals VILMA2's `ρi·ΔI_g` plus the fill (a grounded marine column of
thickness `H` in water depth `d` that melts: `-ρi·H + ρw·d = -ρi·Haf`), up
to the flotation rule. It has no term for the water displaced by bed
motion, `-ρw·C·u`. VILMA2's total load is
`σ = ρi·ΔI_g + ρw·(C·S - F)` (`vilma_sle.f90:310, 433`), so the part that
FastIsostasy represents is

    σ_FI = σ + ρw·C·u = ρi·ΔI_g + ρw·(C·(N + Δφ) - F)

**Channels.** VILMA2's rheology is linear (Maxwell), so its response is
linear in the load history. VILMA2 runs one extra response channel per
FastIsostasy domain: the same Earth structure, its own viscoelastic
memory, driven at each sub-step by `r_d·σ_FI` (`r_d` the weight of `R_d`).
The load is given, so no sea-level iteration is needed. Then
`u_far,d = u - u_d` is the response to everything outside `R_d`, plus the
water displaced by bed motion inside it, which FastIsostasy lacks.

- Linearity is exact for `scheme = "fe"` (the default): each channel is
  advanced, saved and restored in lockstep with every sub-step, including
  the stepper's guard rollbacks (`vilma_timestep.f90:192-227`). Under
  `trap`, the memory endpoint is iterated to a tolerance, so linearity holds
  only to that tolerance. Channels require `fe` (checked at init).
- `u` includes the rotational deformation (`vilma_sle.f90:327`), so
  `u_far` carries the rotational response to the near load. Degree 0 has
  no bed response; the degree-1 part removed with `r_d·σ_FI` (mm to cm) is
  neglected.
- `R_d` is FastIsostasy's active region (`maskactive`). FastIsostasy must
  run with `interactive_sealevel = .true.`, so that its seawater term is in
  `σ_FI`. Checked at init.
- Cost per FastIsostasy domain: one response evaluation and memory advance
  per sub-step (no sea-level iterations), and one more memory state (large
  with 3D viscosity).
- What remains inconsistent: the near field uses FastIsostasy's Earth and
  the far field VILMA2's; the sea surface is VILMA2's everywhere, including
  its coarse geoid response to the near load.

## VILMA2 host API (as is)

`use vilma`, type `solid_earth` (`src/vilma_coupling.f90`):

- `vilma_par_load(se%par, filename, defaults_file, group)`: the `&vilma`
  group overlaid on a complete defaults file.
- `solid_earth_init(se, z_bed_eq, h_ice_eq, grid, h_ice_init)`: global
  reference bed and ice on the host grid. A lon-lat host grid is mapped to
  the Gauss grid by VILMA2 (conservative in, bilinear back).
- `solid_earth_spinup(se, h_ice_lgm)`: optional LGM-memory spin-up.
- `solid_earth_update(se, h_ice, dt_yr)`: advance by `dt_yr`, load ramped
  linearly. Returns `se%rsl`, `se%z_bed = z_bed_eq - rsl`, `se%bsl`.
- `vilma_restart_write/read`, `solid_earth_finalize`.
- Precision `wp = real64`; depends on fesm-utils `coords`, FFTW and SHTns.

Gaps:

1. Only `rsl = N - u + Δφ` reaches the host; `u_lm`/`N_lm` exist inside the
   sea-level solve but are not exported.
2. `se%bsl` is the ice-volume equivalent, not the ocean mean of `C·rsl`
   (`update_bsl`, `vilma_coupling.f90:573`).
3. Densities are hard-coded (931, 1000).
4. Flotation and coastline are decided on the Gauss grid (binary `C`,
   `ocean_function`); at `lmax = 170` that grid is 342 × 680 (about
   0.53°). A cell at a marine margin averages grounded ice with open ocean
   and can be judged afloat, dropping grounded mass from the budget.
5. No response channels.
6. Only executables are built; no `libvilma.a`.
7. README asks for fesm-utils `coords-dev`; yelmox's fesm-utils dev has
   `coords`. Needs a build check.
8. The 3D viscosity reader expects the isostasy_data file; FesmData's
   `Bagge2021_visc.nc` has `log10_visc` on `depth`, 8 structures.

## Parameters

```
[sealevel]
    method         = "vilma2"         ! "const" | "file" | "barystatic" | "vilma2"
    dt             = 10.0             ! [yr] global-layer cadence (as dt_isos)
    grid           = "GLOBAL-1DEG"    ! global grid (FesmData lon-lat)
    bsl_init       = 0.0              ! const
    sl_path        = "none"           ! file
    sl_name        = "none"
    z_bed_ref_file = "fesmdata/Global/{grid}/{grid}_TOPO-<global product>.nc"
    ice_ext        = "file"           ! "ref" | "file": ice outside the footprints
    ice_ext_file   = "fesmdata/Global/{grid}/{grid}_IceHistory-ICE6G_C.nc"
    ice_ext_name   = "H_ice"
    ice_ext_tref   = 0.0              ! [yr] time of the history's present day
    dt_ice_ext     = 100.0            ! [yr] read and interpolation interval
```

`A_ocean_pd` and `A_ocean_path` go away in phase 1b (the ocean area comes
from the global field); phase 1a still needs them.

- `[isostasy<sfx>] method` per domain replaces `[comps] with_isostasy`. The
  FastIsostasy group `[isos<sfx>]` is unchanged.
- `[vilma]`: the full VILMA2 group, read with
  `defaults_file = "input/vilma_defaults.nml"`; `scripts/check_par_nml`
  validates it. The par files set the CLIMBER-X overrides (`lmax = 170`,
  the PREM radial structure, `l_visc_3d = .true.`), and
  `visc_3d_file` / `visc_3d_structure` point to FesmData's
  `Bagge2021_visc.nc`, structure 6. The densities are set by yelmox, not
  in `[vilma]`.
- `[sealevel] grid` sets the coastline resolution outside host regions;
  VILMA2's pointwise products are on its Gauss grid (about 0.53° at
  `lmax = 170`). See open question 3.

## Restart and output

- Run bundle: `sealevel_restart.nc` (`barystatic`: `Δφ`, bsl, `C`) or
  `vilma_restart.nc` (`vilma2`: VILMA2's state, the channel memories and
  the last host terms). `const` and `file` have no state. `bsl_restart.nc`
  is retired in phase 1b.
- Domain bundle: `isos_restart.nc` under `fastisostasy`. Restart files from
  before phase 1b carry the unassembled reference ice and fail the
  reference check, so phase 1b needs cold starts.
- Run level: `sealevel.nc` (global grid: ice, `I_g`, `C`, host weights;
  under `vilma2` also `dz_bed`, `dz_ss`, `rsl`) and `sealevel_ts.nc`
  (`bsl`, `bsl_ive`, `Δφ`, ocean area, contribution per domain and
  external, clipped mass; under `vilma2` the worst mass residual and solver
  time).
- Domain: `isos.nc` / `isos_ts.nc` with the neutral fields (`z_bed`,
  `z_ss`, `dz_bed`, `dz_ss`, `bsl`); `w`, `we`, `u_far` under
  `fastisostasy`.

## Build

- `vilma=1` (default): `VILMAROOT` (configme dependency),
  `INC_VILMA`/`LIB_VILMA`, plus SHTns from fesm-utils (`make shtns`); FFTW
  is already needed by FastIsostasy. Compiles `libs/sealevel_vilma.f90`,
  the adapter that owns the VILMA2 state.
- `vilma=0`: `libs/sealevel_vilma_stub.f90` with the same interface;
  `method = "vilma2"` stops at init with "rebuild with vilma=1".
- `libvilma.a` is a file target, and host objects depend on it.

## Changes outside yelmox

VILMA2:

1. Export `dz_bed = u` and `dz_ss = N + Δφ` on the host and Gauss grids.
2. `se%bsl` = ocean mean of `C·rsl`; the current diagnostic is renamed
   (ice-volume equivalent).
3. Densities as parameters; the GIA benchmark values stay VILMA2's
   defaults.
4. Host-resolved terms: optional host-grid inputs `w`, `I_g`, `C`, `F` at
   init (reference) and update. They are ramped linearly over the update
   interval like the ice (the previous values are the ramp start and go
   into the restart; the `dt = 0` seed step uses the current values).
   Inside `w` the sea-level solve uses them; on `1 - w` it applies its own
   flotation and coastline rule, and the outer coastline loop iterates only
   that part. `C` and `w` are not rescaled by the mass-conserving remap.
   With `w ≡ 0` results are bit-identical to now.
5. Response channels: region weights registered at init; each channel
   driven by `r_k·σ_FI`, advanced, saved and restored in lockstep with the
   stepper, in spin-up too; `u_k` exported; memories in the restart;
   `scheme = "fe"` required.
6. Reader for FesmData `Bagge2021_visc.nc` (`log10_visc`, `depth`,
   `visc_3d_structure`).
7. `make vilma-static` → `libvilma/include/libvilma.a` + `.mod`, `openmp`
   variants.
8. Build and tests against fesm-utils dev; README.

FastIsostasy:

1. `isos_update` and `isos_init_state` take the uniform sea-surface shift
   as a scalar and no longer write bsl back. The restart gate stays inside
   FastIsostasy: at the redo the passed shift is ignored and the restored
   `z_ss` kept (`fastisostasy.f90:829-834, 849`).
2. Phase 1a only: the sea-level contribution (`calc_sl_contributions`) as
   a separate call before `isos_update`, also at init, with the same gate
   (`use_restart` and `time <= time_prognostics`, `fastisostasy.f90:801`).
3. Optional external fields: a sea-surface anomaly (replaces its regional
   geoid) and a bed displacement (enters flotation and the ocean column,
   not the solve state).

FesmData:

1. A global present-day topography product (`z_bed`, `H_ice`, `z_srf`) on
   the Global lon-lat grids: the domain merges (BedMachine Antarctica v4
   and Greenland v6, IceBoost v2.0, GEBCO 2026) applied globally.
2. `IceHistory-ICE6G_C` on the Global grids (exists).

Yelmo: none. kryos calls Yelmo's flotation routines for `f_grnd_ref` and
for the external ice in host regions.

## Phases

- **0 (VILMA2, FesmData):** VILMA2 changes 1-8 and the FesmData product.
  Tests: benchmarks still pass; `rsl = dz_ss - dz_bed` to round-off;
  `w ≡ 0` bit-identical; host terms from a fine grid converge to VILMA2's
  own flotation on the same geometry (same densities) as the Gauss grid is
  refined; under `fe`, channels plus remainder equal the total response to
  round-off; the FesmData 3D file matches the isostasy_data file for
  structure 6.
- **1a (yelmox, FastIsostasy):** `kryos_sealevel` and `kryos_isostasy` with
  `const`, `file`, and `barystatic` provisionally computed as today:
  FastIsostasy's volume contributions (init included) summed by the run on
  each domain's `dt_isos`, `A_ocean_now` from the bsl before the
  contribution, `Δφ ≡ bsl`. FastIsostasy changes 1-2. Par files get an
  explicit `[sealevel]` group with their current behaviour (`fastiso` →
  `barystatic`). Bit-identical on the single-domain Levante cases; bipolar
  changes only by the order fix.
- **1b (yelmox):** the global ice field (assembly, footprints, host-resolved
  terms, reference from `dta%pd`), the `barystatic` budget replaces the
  provisional one, FastIsostasy's load comes from the assembly,
  `sealevel.nc`, `sealevel_restart.nc`. Validated against 1a: bsl and ice
  volume on the ANT, GRL and bipolar cases. Expected differences: the
  budget formula against the volume accounting, densities, `∫C` against
  the `A_ocean(bsl)` table, spherical against projected areas, the `Δφ`
  correction, the masked-reference fix, the one-step lag of `u_R`.
- **2 (yelmox):** VILMA2 adapter and stub, `vilma=1`,
  `input/vilma_defaults.nml`, `[vilma]`, `[isostasy] method = "global"`,
  restart, output.
- **3 (yelmox, FastIsostasy):** FastIsostasy with VILMA2: FastIsostasy
  change 3, response channels, `u_far`.
- **4 (tests):** VILMA2 standalone vs yelmox with prescribed ice (same
  load); ANT deglaciation under `fastisostasy`, `global` and both;
  bipolar (two domains, one Earth, two channels); exact restart for every
  method; cost of the default (T170, 3D viscosity, channels).

## Open questions

1. Footprint: `MASK_ICE_DYNAMIC` only, or every cell with
   `mask_ice /= MASK_ICE_NONE` (also `FIXED`)? With dynamic only, fixed
   cells take the external ice rather than the ice Yelmo prescribes there
   (`ice_bnd`), so Yelmo's geometry and the load differ wherever `ice_bnd`
   and `ice_ext` use different histories.
2. Default ice history: ICE-6G_C covers 26-0 ka only, so runs that start
   earlier stop at init. Longer ones in FesmData: PaleoMIST 1.0 (global,
   80-0 ka), GLAC-1D nn9927 (120-0 ka, no Eurasia). Or default to `"ref"`.
3. Default `[sealevel] grid`: `GLOBAL-1DEG`, or `GLOBAL-0.5DEG` to match
   VILMA2's Gauss grid at `lmax = 170`?
4. When do the existing par files switch from `barystatic` to `vilma2`:
   with phase 2, or after the phase-4 validation?
5. Reference equilibrium: the present-day bed treated as relaxed (as
   FastIsostasy does; ongoing adjustment only from the spin-up memory), or
   a relaxed ice-free bed as in CLIMBER-X (`z_bed_rel_file`)?
