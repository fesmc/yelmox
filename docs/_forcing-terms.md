# Forcing as terms: series × pattern

Status: design (2026-10-21), agreed on the interactive mock-up
`docs/_forcing-terms-mockup.html` (open it in a browser; the series are the
real files of `input/`). Nothing implemented. The leading `_` keeps quarto
from rendering it. Supersedes sections 3, 4, 6, 7 and 8 of
`_forcing-design.md` (indices, axes, modifiers, OBM); the rest of that doc
(insolation, database, retirement of snapclim, phases 0-6 done) stands.

State checked: yelmox dev 54d9ee6 (`libs/snapshots.f90`, `libs/climate_mods.f90`,
`libs/kryos_forcing.f90`, `input/forcing_db.nml`, `input/yelmox_defaults.nml`).

## 1. Why

The concepts of `_forcing-design.md` are sound; their encoding is not. One
physical statement, "blend piControl → LGM by the glacial index, temperature
additively, precipitation as a ratio", is spread over `[indices] names`,
`[idx_a]`, `[snap] snapshots`, `[snap] axes`, `[snap_axis_a]`,
`[snap_field_tas]`, `[snap_field_pr]`: four name-list + per-name-group
registries with the same `snap_` prefix, and the uniform anomaly `[anom]`,
the ocean rules `[ocn_fraction]` and `[ocn_obm]` as unrelated groups with
their own vocabulary. Names such as `a`, `at`, `one` carry assumptions
(orbital, temperature) that the equation does not.

The fix: write the forcing as one equation, and make the par file a
transcription of it. The names are mathematical (a series α(τ), a term, a
pattern), the user decides what α is.

## 2. The equations

The forcing at one moment, snapshot backend. x cell, m month, z depth,
t model time.

```
0. Forcing time
   τ(t)  = time_rel, set by [ctrl] tstep_method:
           const: τ = tstep_const      rel: τ = t      cal: calendar(t)

1. Scalar series  α_k(τ), named by the user (alpha, beta, gmt, dTo, ...)
   file         a table or netcdf series, optionally normalised:
                α = (F(τ) − F̄[t0]) / (F̄[t1] − F̄[t0]),  t0, t1 a time or a window (mean)
   const        α = f_min
   generated    ramp-slope, ramp-time, sin, ... (tsgen)
   feedback     PI42, H312PID, ... on the response V_ice (tsgen)
   derived      computed each step from the model state (section 3.2):
                the mean sea-level air anomaly, an OBM box temperature

2. Snapshots S_j (named sets of fields, from the library), reduced to sea level once
   T_j^sl(x,m) = T_j + Γ_m z_j         P_j^sl(x,m) = P_j · exp(β_p Γ_m z_j)
   Θ_j(x,z), S_j(x,z) as given.        S_j may depend on τ (a reconstruction in time).
   S_0 is the reference.

3. Composition at sea level: reference + Σ_k series × pattern
   T^sl(x,m,τ) = T_0^sl + Σ_k α_k(τ) [T_{k1}^sl − T_{k0}^sl]  + Σ_u f_ta,u α_u(τ)
   P^sl(x,m,τ) = P_0^sl · Π_k [1 + α_k(τ) (P_{k1}^sl / P_{k0}^sl − 1)] · Π_u exp(β_p f_ta,u α_u(τ))
   Θ(x,z,τ)    = Θ_0 (or the series of a replace term) + Σ_k α_k [Θ_{k1} − Θ_{k0}] + Σ_u f_to,u α_u
   S(x,z,τ)    = S_0 (or ...)                          + Σ_k α_k [S_{k1} − S_{k0}] + Σ_u f_so,u α_u
   k: terms with a snapshot pattern; u: terms with a uniform pattern.
   A term without snap0 adds snap1 itself (a given anomaly or ratio).
   Order: the atmosphere, then the derived series, then the ocean.

4. To the model surface z_s(x,t)
   T(x,m) = T^sl − Γ_m z_s            P(x,m) = P^sl · exp(−β_p Γ_m z_s)

5. Static transforms (no series)
   taper:  X ← X_0(z_s) + w(x) (X − X_0(z_s)) [+ (1−w) ⟨X − X_0⟩_{w=1}],  X ∈ {T, P}
   bias:   T ← T + f_tas b_T(x,m[,τ]),   P ← P (1 + f_pr (b_P(x,m[,τ]) − 1))

6. Also at the forcing time
   insolation  I(φ, d; orbit(τ))  |  orbit(t_fixed)  |  orbit(e, ε, ϖ)
   sea level   η(τ) from a file  |  η_0  |  from the ice volume
   border ice  H_b(x) = H(x, τ)  |  α_k(τ) H(x, τ_slice)

Consumers unchanged: surface (T, P, I → smb, T_srf), shelf (Θ(z_base), S → melt).
```

Rules fixed by field, not configurable: `tas`, `to`, `so` compose additively;
`pr` as a ratio, and in a uniform term follows the term's `f_ta` as
exp(β_p ΔT). Γ (annual, summer; a monthly cosine between them) and β_p are
`[snapshots] lapse` and `f_p`, the same values in the sea-level reduction,
the uniform term and the elevation transform. A user who wants no
precipitation response to a uniform term leaves `pr` out of its fields.

Backend variants. **esm** replaces 2-3 by its own calendar-time anomaly and
refuses terms on `tas`/`pr` when it supplies the smb (`surface_method =
"climate"`); ocean terms apply at the shelf base. **rembo** takes the uniform
`tas` term as its input: ΔT_sum = f_ta α, ΔT_ann = 1.3 ΔT_sum; its ocean from
the snapshot backend as now.

## 3. Concepts and contracts

### 3.1 Series

A named scalar on the forcing time, owned by the driver (`indices_class`
today), advanced once per step, restarted from `series_<name>_restart.nc`,
written to the 1D output. Shared by every domain of a multi-domain driver.
`with_kill` and `restart_every_df` per series as now. The same series may
drive several terms.

### 3.2 Derived series

A series computed each step by the climate step, not advanced by tsgen:

| `from` | value | keys |
|---|---|---|
| `air_anomaly` | factor × mean over `region` of (T^sl_ann − T^sl_ann,0) of this update | `factor`, `region` |
| `obm_north`, `obm_south` | the box temperature θ, or θ − θ_init with `anomaly = True` | `anomaly` |

Contract: a derived series is evaluated after the atmosphere is composed,
so a term driven by it may touch `to` and `so` only (an error at init
otherwise). `air_anomaly` needs a backend with `tsl_ann` (snapshots); the
OBM series exist in `yelmox_bipolar` only, their value set by the bridge
where `obm2ism` runs today. No lag: the anomaly of the same update (as
decided in `_forcing-design.md` section 7). A derived series is not
restarted (recomputed) and is written to the 1D output like any series.

### 3.3 Snapshots

As today: named groups binding fields to varslice groups, static or in time,
reduced to sea level at loading; the library `input/forcing_db.nml`
(`[snapshot_<name>]`), the fields in `input/clim_db.nml`. A par file may add
or override a snapshot. The reference is `[forcing] reference`. The set of
snapshots loaded is derived from the reference and the terms; no list.

### 3.4 Terms

One term = one series × one pattern, on an explicit list of fields:

| `pattern` | adds | keys |
|---|---|---|
| `snapshots` (default) | α (snap1 − snap0), or α snap1 when `snap0 = ""` | `snap0`, `snap1`, `fields` |
| `uniform` | f_ta α to `tas` (and exp(β_p f_ta α) to `pr`), f_to α to `to`, f_so α to `so` | `fields`, `f_ta`, `f_to`, `f_so` (default 1) |
| `replace` | the field = α at every depth, in place of the reference; other terms add to it | `fields` (`to`, `so` only), one per field |

`series = ""` means the constant 1 (the full anomaly of a reconstruction).
Contracts, checked at init with an error naming the term: every field is
bound by `snap0` (if given), `snap1` and the reference; `snap1` of a term
without `snap0` is a sea-level anomaly or ratio (no `zs`); `replace` once
per field; a derived series on the ocean only.

### 3.5 Transforms and other consumers

Elevation (always), `[taper]` and `[bias]` as today; the bias groups may
carry a time axis and are then sliced at τ (varslice, as time-varying
snapshots). `[orbit] method = "forcing_time" | "fixed_time" | "fixed_orbit"`
(today `clock | time | fixed`; with `fixed_orbit` the three elements are
shown). Sea level and `[ice_bnd]` unchanged. In the docs the "forcing clock"
becomes the "forcing time" τ (`time_rel`), set by `[ctrl] tstep_method`,
which stays as it is (it defines what time means; the forcing then uses it).

## 4. The par file

One database list, `[forcing] var_defs`, searched by every consumer
(series, snapshots and their fields, taper, bias, ice_bnd): the par file
first, then the files in order (the lookup rule of `_par-files-plan.md`).
Unset keys from `input/yelmox_defaults.nml` (templates `series`, `term`,
`snapshot`).

```fortran
&snapshots
    lapse = 8.0e-3 6.5e-3     ! [K/m] Γ annual, summer
    f_p   = 0.06              ! [1/K] β_p
/

&forcing
    reference = "pd"
    series    = "alpha" "dTo"
    terms     = "glacial" "ocean"
    var_defs  = "input/forcing_db.nml" "input/clim_db.nml"
/

&series_alpha
    method      = "series"
    series_file = "input/alpha_combined_125kyr_interp.dat"
/

&series_dTo
    method = "derived"
    from   = "air_anomaly"
    factor = 0.25
    region = "all"
    units  = "K"
/

&term_glacial
    series = "alpha"
    snap0  = "picontrol"
    snap1  = "lgm"
    fields = "tas" "pr"
/

&term_ocean
    series  = "dTo"
    pattern = "uniform"
    fields  = "to"
/
```

GRL paleo A reads: T = T_pd + α (T_lgm − T_pic), P = P_pd [1 + α (P_lgm/P_pic − 1)],
Θ = Θ_pd + 0.25 ⟨T^sl − T^sl_pd⟩. The other cases (the mock-up presets):

| Case | series | terms |
|---|---|---|
| GRL paleo B | `dTo` | `buizert` (1, b18_base → b18(τ), tas); `badgeley` (1, "" → b20(τ), pr); `ocean` |
| LIS, Banderas 2018 | `alpha` const 0.87, `beta` const 0.48 (or NGRIP), `dTo` | `orbital` (α, picontrol → lgm), `millennial` (β, lgmweak → lgmstrong), `ocean` |
| ANT paleo | as GRL A on `pd` (ANT) | |
| Pyrenees | `gamma` orbital, normalised 0 at 0, 1 at −21 ka | `cooling` (γ, uniform, tas pr, f_ta = −10) |
| Hysteresis | `gamma` ramp or PI42 | `warming` (γ, uniform, tas pr to, f_to = 0.25) |
| Bipolar north | `alpha`, `gmt` [K], `tn` derived obm_north | `glacial` (α), `regional` (gmt, uniform, f_ta = 1.07), `ocean` (tn, replace, to) |
| Present day | none | none |

A term is one group in the library or the par file; an experiment adds or
removes terms with `-p forcing.terms==glacial,ocean`, which is where the
par-files plan's `-p` arrays shrink to.

## 5. From today's groups

| today | new |
|---|---|
| `[indices] names`, `var_defs` | `[forcing] series`, `[forcing] var_defs` (the one list) |
| `[idx_<name>]` | `[series_<name>]` (+ `method = "derived"`, `from`, `factor`, `region`, `anomaly`) |
| `[snap]` | `[snapshots]`: `lapse`, `f_p`, `grid_global`, `combine`; `ref_name` → `[forcing] reference`; `snapshots`, `axes`, `fields` dropped (derived from the terms) |
| `[snap_<name>]` | `[snapshot_<name>]` |
| `[snap_axis_<name>]` `lo`, `hi`, `index` | `[term_<name>]` `snap0`, `snap1`, `series`, `fields` |
| `[snap_field_<name>]` `blend`, `axes` | dropped: the rule is fixed by field, the fields are listed by the term |
| `[anom]` | a `uniform` term; rembo reads the uniform `tas` term |
| `[ocn_fraction]` | a derived series `air_anomaly` + a `uniform` term on `to` |
| `[ocn_obm]` absolute / anomaly | a derived series `obm_<hemisphere>` (`anomaly`) + a `replace` / `uniform` term on `to`, applied by the bridge |
| `[obm_forcing]` | unchanged (a consumer of a series, as `[ice_bnd]`) |
| `[taper]`, `[bias]`, `[orbit]`, `[barysealevel]`, `[ice_bnd]` | unchanged; `[bias]` groups may have a time axis; `[orbit] method` renamed |
| `idx_one` | gone (`series = ""`) |
| bipolar `snap_north_pd`, `snap_north_axis_at`, ... | `[forcing_north]` selecting library snapshots by content name (`pd_grl`, `pd_ant`); `[snapshots_north]` |

Library names become descriptive and free of use: `series_glacial_alpha`
rather than `idx_a`; the par file names what it is in this run.

## 6. Design idea, not planned: the equation as text

```fortran
&forcing
    series = "alpha" "beta"
    tas    = "pd + alpha*(lgm - picontrol) + beta*(lgmstrong - lgmweak)"
    pr     = "pd * (1 + alpha*(lgm/picontrol - 1)) * (1 + beta*(lgmstrong/lgmweak - 1))"
    to     = "pd + dTo"
/
```

The most readable form, and what the mock-up shows as tab (b). Not pursued:
a parser in Fortran, and the precipitation line either repeats what the
temperature line implies or must be derived by a per-field rule, which is
what the terms do explicitly. If wanted later, a reader for this grammar
can produce the same term structure.

## 7. Code

| File | Change |
|---|---|
| `libs/kryos_forcing.f90` | `indices_class` → `series_class`: `[forcing] series`, `[series_<name>]`; kind `derived` (no tsgen; value set by `series_set(name, value)`), restart and output as now |
| `libs/snapshots.f90` | `[snapshots]`, `[forcing] reference`, terms read from `[term_<name>]` (snapshot pattern); the snapshot set derived from the terms; `combine` per field by the fixed rule; `snap_field` registry → the term's field list |
| `libs/climate_mods.f90` | `anom`, `ocn_fraction`, `ocn_obm` removed; `terms_apply_uniform` (uniform and replace terms on `climate_out`), `taper`, `bias` (time axis) stay; `derived_update` computes `air_anomaly` after the atmosphere |
| `libs/yelmox_climate.f90` | step order: backend (snapshot terms) → uniform `tas`/`pr` terms → derived series → ocean terms → taper → bias; contract checks at init |
| `yelmox_bipolar/obm_coupling.f90` | sets the `obm_*` derived series from the boxes where `obm2ism` ran; `[ocn_obm_<sfx>]` removed |
| `libs/climate_rembo.f90` | reads the uniform `tas` term instead of `[anom]` |
| `input/yelmox_defaults.nml` | templates `series`, `term`, `snapshot`, `snapshots`, `forcing`; the old groups removed |
| `input/forcing_db.nml` | `snapshot_*`, `series_*`, `term_*` with descriptive names |
| `yelmox/*.nml`, `yelmox_bipolar/*.nml`, `scripts/*/params.sh` | converted; the `-p` arrays shrink to `forcing.terms` plus run control |
| `scripts/check_par_nml.py` | templates by prefix `series_`, `term_`, `snapshot_` |
| `scripts/forcing-builder/`, `docs/forcing-builder/` | export the new groups; the equation panel of the mock-up moves into the builder |
| docs | `yelmox.md`, `climate-snapshots.md`, experiment pages, CHANGELOG; "forcing time" wording |

## 8. Phases

All on dev, one branch `forcing-terms`, merged when the whole chain passes.

| # | Content | Check |
|---|---|---|
| 1 | this doc, the mock-up in `docs/` | |
| 2 | series registry: rename, `derived` kind, `series_set`; configs unchanged through an alias of the old keys is **not** done: the configs convert in 4 | unit: a derived series set and read |
| 3 | terms in the snapshot backend, uniform and replace terms, derived `air_anomaly`, bridge sets the OBM series; `anom`/`ocn_fraction`/`ocn_obm` removed | |
| 4 | library, defaults, bases, bipolar, scripts, checker | `check_par_nml` clean |
| 5 | docs-figures 18 cases + the bipolar case on Levante | bit-identical to dev 54d9ee6 (the same operations in the same order) |
| 6 | forcing builder on the new groups, with the equation panel; docs; CHANGELOG | |

## 9. Decisions taken (2026-10-21)

Encoding (a), one group per term; (b) recorded in section 6. The uniform
anomaly is a term. Derived series as a series kind, with the factor in the
series. Fields always explicit. Precipitation by the fixed rule, nothing to
set. `snap0`/`snap1` for the two snapshots of a pattern (0 and 1 as the
index values, tied to `norm_t0`/`norm_t1`); `f_ta/f_to/f_so` kept;
`tstep_method` kept; `forcing_time` wording, `clock` dropped. OBM absolute
mode as the `replace` pattern. Time-dependent bias through the varslice time
axis, not as a term. Normalisation windows kept (`norm_t0 = t` or `t_a t_b`).
Taper, bias, orbit, sea level, border ice stay outside the sum.

Open: the group names `[forcing]`, `[snapshots]`, prefixes `series_`,
`term_`, `snapshot_` (as on the mock-up unless changed); whether
`[snapshots] combine` keeps a use.
