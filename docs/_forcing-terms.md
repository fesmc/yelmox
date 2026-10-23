# Forcing as terms: coefficient × pattern

Status: design agreed 2026-10-21 on the interactive mock-up
`docs/_forcing-terms-mockup.html` (open it in a browser; the records are the
real files of `input/`); phases 2-4 implemented on the branch `forcing-terms`
(2026-10-21, not yet verified on Levante). The leading `_` keeps quarto
from rendering it. Supersedes sections 3, 4, 6, 7 and 8 of
`_forcing-design.md` (indices, axes, modifiers, OBM); the rest of that doc
(insolation, database, retirement of snapclim, phases 0-6 done) stands.
2026-10-22: series → coefficient (`[forcing] coefs`, `[coef_<name>]`, `[term] coef`), the vocabulary coefficient / pattern / parameter / term agreed on the forcing builder.

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
transcription of it. The names are mathematical (a coefficient α(τ), a
term, a pattern, a parameter), the user decides what α is.

## 2. The equations

The forcing at one moment, snapshot backend. x cell, m month, z depth,
t model time.

```
0. Forcing time
   τ(t)  = time_rel, set by [ctrl] tstep_method:
           const: τ = tstep_const      rel: τ = t      cal: calendar(t)

1. Coefficients  α_k(τ), scalars named by the user (alpha, beta, gmt, dTo, ...)
   series       a table or netcdf time series, optionally normalised:
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

3. Composition at sea level: reference + Σ_k coefficient × pattern
   T^sl(x,m,τ) = T_0^sl + Σ_k α_k(τ) [T_{k1}^sl − T_{k0}^sl]  + Σ_u f_ta,u α_u(τ)
   P^sl(x,m,τ) = P_0^sl · Π_k [1 + α_k(τ) (P_{k1}^sl / P_{k0}^sl − 1)] · Π_u exp(β_p f_ta,u α_u(τ))
   Θ(x,z,τ)    = Θ_0 (or the coefficient of a replace term) + Σ_k α_k [Θ_{k1} − Θ_{k0}] + Σ_u f_to,u α_u
   S(x,z,τ)    = S_0 (or ...)                          + Σ_k α_k [S_{k1} − S_{k0}] + Σ_u f_so,u α_u
   k: terms with a snapshot pattern; u: terms with a uniform pattern.
   A term without snap0 adds snap1 itself (a given anomaly or ratio).
   Order: the atmosphere, then the derived coefficients, then the ocean.

4. To the model surface z_s(x,t)
   T(x,m) = T^sl − Γ_m z_s            P(x,m) = P^sl · exp(−β_p Γ_m z_s)

5. Static transforms (no coefficient)
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
the uniform term and the elevation transform. These and the gains `f_ta`,
`f_to`, `f_so` of a uniform term are the parameters: fixed numbers set in the
par file. A user who wants no
precipitation response to a uniform term leaves `pr` out of its fields.

Backend variants. **esm** replaces 2-3 by its own calendar-time anomaly and
refuses terms on `tas`/`pr` when it supplies the smb (`surface_method =
"climate"`); ocean terms apply at the shelf base. **rembo** ignores `pr` in a
uniform term (a notice at init) and takes the uniform
`tas` term as its input: ΔT_sum = f_ta α, and REMBO's annual mean
ΔT_ann = (1 + T_wintfac)/2 ΔT_sum (its seasonal cycle, a diagnostic); its ocean from
the snapshot backend as now.

## 3. Concepts and contracts

### 3.1 Coefficients

A named scalar on the forcing time, owned by the driver (`coefs_class`),
advanced once per step, restarted from `coef_<name>_restart.nc` (an older
`series_<name>_restart.nc` is still read, with a note), written to the 1D
output as `coef_<name>`. Shared by every domain of a multi-domain driver.
`with_kill` and `restart_every_df` per coefficient as now. The same
coefficient may drive several terms. Its kind is the tsgen `method`: `series`
(a time series from a file), `const`, a ramp or wave, a controller, or
`derived`.

### 3.2 Derived coefficients

A coefficient computed each step by the climate step, not advanced by tsgen:

| `from` | value | keys |
|---|---|---|
| `air_anomaly` | factor × mean over `region` of (T^sl_ann − T^sl_ann,0) of this update | `factor`, `region` |
| `obm_north`, `obm_south` | the box temperature θ, or θ − θ_init with `anomaly = True` | `anomaly` |

Contract: a derived coefficient is evaluated after the atmosphere is composed,
so a term driven by it may touch `to` and `so` only (an error at init
otherwise). `air_anomaly` needs a backend with `tsl_ann` (snapshots); the
OBM coefficients exist in `yelmox_bipolar` only, their value set by the bridge
where `obm2ism` runs today. No lag: the anomaly of the same update (as
decided in `_forcing-design.md` section 7). A derived coefficient is not
restarted (recomputed) and is written to the 1D output like any coefficient.

### 3.3 Snapshots

As today: named groups binding fields to varslice groups, static or in time,
reduced to sea level at loading; the library `input/forcing_db.nml`
(`[snapshot_<name>]`), the fields in `input/clim_db.nml`. A par file may add
or override a snapshot. The reference is `[forcing] reference`. The set of
snapshots loaded is derived from the reference and the terms; no list.

### 3.4 Terms

One term = one coefficient × one pattern, on an explicit list of fields:

| `pattern` | adds | keys |
|---|---|---|
| `snapshots` (default) | α (snap1 − snap0), or α snap1 when `snap0 = ""` | `snap0`, `snap1`, `fields` |
| `uniform` | f_ta α to `tas` (and exp(β_p f_ta α) to `pr`), f_to α to `to`, f_so α to `so` | `fields`, `f_ta`, `f_to`, `f_so` (default 1) |
| `replace` | the field = α at every depth, in place of the reference; other terms add to it | `fields` (`to`, `so` only), one per field |

`coef = ""` means the constant 1 (the full anomaly of a reconstruction).
Contracts, checked at init with an error naming the term: every field is
bound by `snap0` (if given), `snap1` and the reference; `snap1` of a term
without `snap0` is a sea-level anomaly or ratio (no `zs`); `replace` once
per field, on ocean profiles; a uniform term on `to`/`so` needs an ocean
(profiles or at the shelf base); a derived coefficient on the ocean only.
One predicate of the backend serves the uniform terms on `tas`/`pr` and the
transforms (`[taper]`, `[bias]`), which change the atmosphere after the
backend: its atmosphere must drive the surface (monthly `tas` and `pr`, now
and reference). They are refused when the backend supplies the smb itself
(esm with its smb; rembo, which takes the uniform `tas` terms as its input).

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
(coefficients, terms, snapshots and their fields, taper, bias, ice_bnd): the par
file first, then the files in order (the lookup rule of
`_par-files-plan.md`); default `"input/clim_db.nml" "input/forcing_db.nml"`.
Unset keys from `input/yelmox_defaults.nml` (templates `forcing`, `coef`,
`term`, `snapshots`, `snapshot`).

The coefficients are driver-owned and shared by every domain, the reference and
the terms belong to a domain: in `yelmox` all four keys sit in one
`[forcing]`; in `yelmox_bipolar` `[forcing]` holds `coefs` and `var_defs`,
`[forcing_north]` and `[forcing_south]` hold `reference` and `terms`, and
the snapshots are named by content (`snapshot_pd_grl`, `snapshot_pd_ant`)
since the library has no domain suffix.

```fortran
&snapshots
    lapse = 8.0e-3 6.5e-3     ! [K/m] Γ annual, summer
    f_p   = 0.06              ! [1/K] β_p
/

&forcing
    reference = "pd"
    coefs     = "alpha" "dTo"
    terms     = "lgm" "ocean"
    var_defs  = "input/forcing_db.nml" "input/clim_db.nml"
/

&coef_alpha
    method      = "series"
    series_file = "input/alpha_combined_125kyr_interp.dat"
/

&coef_dTo
    method = "derived"
    from   = "air_anomaly"
    factor = 0.25
    region = "all"
    units  = "K"
/

&term_lgm
    coef   = "alpha"
    snap0  = "picontrol"
    snap1  = "lgm"
    fields = "tas" "pr"
/

&term_ocean
    coef    = "dTo"
    pattern = "uniform"
    fields  = "to"
/
```

GRL paleo A reads: T = T_pd + α (T_lgm − T_pic), P = P_pd [1 + α (P_lgm/P_pic − 1)],
Θ = Θ_pd + 0.25 ⟨T^sl − T^sl_pd⟩. The other cases (the mock-up presets):

| Case | coefficients | terms |
|---|---|---|
| GRL paleo B | `dTo` | `buizert` (1, b18_base → b18(τ), tas); `badgeley` (1, "" → b20(τ), pr); `ocean` |
| LIS, Banderas 2018 | `alpha` const 0.87, `beta` const 0.48 (or NGRIP), `dTo` | `lgm` (α, picontrol → lgm), `millennial` (β, lgmweak → lgmstrong), `ocean` |
| ANT paleo | as GRL A on `pd` (ANT) | |
| Pyrenees | `alpha` (t2018 orbital record, normalised 0 at 0, 1 at −21 ka) | `uniform` (α, uniform, tas pr, f_ta = −10) |
| Hysteresis | `gamma` ramp or PI42 | `warming` (γ, uniform, tas pr to, f_to = 0.25) |
| Bipolar north | `alpha` (r2020 orbital record, normalised), `gmt` [K], `tn` derived obm_north | `lgm` (α), `regional` (gmt, uniform, f_ta = 1.07), `ocean` (tn, replace, to) |
| Present day | none | none |

A term is one group in the library or the par file; an experiment adds or
removes terms with `-p forcing.terms==lgm,ocean`, which is where the
par-files plan's `-p` arrays shrink to.

## 5. From today's groups

| today | new |
|---|---|
| `[indices] names`, `var_defs` | `[forcing] coefs`, `[forcing] var_defs` (the one list) |
| `[idx_<name>]` | `[coef_<name>]` (+ `method = "derived"`, `from`, `factor`, `region`, `anomaly`) |
| `[snap]` | `[snapshots]`: `lapse`, `f_p`, `grid_global`; `combine` and `f_p_ne` removed (no effect); `ref_name` → `[forcing] reference`; `snapshots`, `axes`, `fields` dropped (derived from the terms) |
| `[snap_<name>]` | `[snapshot_<name>]` |
| `[snap_axis_<name>]` `lo`, `hi`, `index` | `[term_<name>]` `snap0`, `snap1`, `coef`, `fields` |
| `[snap_field_<name>]` `blend`, `axes` | dropped: the rule is fixed by field, the fields are listed by the term |
| `[anom]` | a `uniform` term; rembo takes the uniform `tas` terms as ΔT_sum (its ocean anomaly `1.3 f_ta f_to α` is gone: a uniform term on `to` instead) |
| `[ocn_fraction]` | a derived coefficient `air_anomaly` + a `uniform` term on `to` |
| `[ocn_obm]` absolute / anomaly | the library coefficients `obm_north`/`obm_south` (derived, `anomaly = False`) + the `replace` terms `term_obm_north`/`_south` (`uniform` with `anomaly = True` for the anomaly mode); the bridge sets the coefficients and re-applies the ocean terms (idempotent) |
| `[obm_forcing]` | `index` → `coef` (a consumer of a coefficient, as `[ice_bnd]`) |
| `[taper]`, `[bias]`, `[orbit]`, `[barysealevel]`, `[ice_bnd]` | unchanged; `[bias]` groups may have a time axis; `[orbit] method` renamed; `[ice_bnd] index` → `coef` |
| `idx_one` | gone (`coef = ""`) |
| bipolar `snap_north_pd`, `snap_north_axis_at`, ... | `[forcing_north]` selecting library snapshots by content name (`pd_grl`, `pd_ant`); `[snapshots_north]` |

Coefficients are named by their symbol in the equation (`alpha`, `beta`,
`gamma`, ...), never by an interpretation of the record; which record a
symbol loads is the setup's, in its par file: `[coef_alpha]` is `alpha_combined_125kyr`
(as is) in the paleo bases, `alpha_orbital_r2020` (normalised 0 at 0, 1 at
21 ka) in `yelmox_Greenland.nml` and the bipolar north, `alpha_orbital_t2018`
(normalised) in `yelmox_Pyrenees.nml`, the constant 1 in LIS and North. The
library keeps a default `coef_alpha` (combined) for a par file without one,
and the physical or shaped coefficients `coef_dTo` (derived, 0.25),
`coef_ramp`, `coef_hyst`, `coef_dT`, `coef_obm_north/south`. Terms
are named by what they add: `term_lgm` (lgm − picontrol on α), `term_b18`,
`term_b20`, `term_ocean`, `term_ramp`, `term_hyst`, `term_dT`,
`term_obm_north/south`. A par file redefines a library group to change it for
its run.

The `-p` arrays of the experiments shrink to `forcing.coefs==...` and
`forcing.terms==...` plus the transforms (scripts/grl-paleo/params.sh,
scripts/ant-paleo/params.sh).

## 6. Design idea, not planned: the equation as text

```fortran
&forcing
    coefs  = "alpha" "beta"
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
| `libs/kryos_forcing.f90` | `indices_class` → `coefs_class`: `[forcing] coefs`, `[coef_<name>]`; kind `derived` (no tsgen; value set by `coefs_set(coefs, name, value)`), restart and output as now |
| `libs/snapshots.f90` | `[snapshots]`, `[forcing] reference`, terms read from `[term_<name>]` (snapshot pattern); the snapshot set derived from the terms; the blend per field by the fixed rule; `snap_field` registry → the term's field list |
| `libs/climate_mods.f90` | `anom`, `ocn_fraction`, `ocn_obm` removed; `terms_apply_uniform` (uniform and replace terms on `climate_out`), `taper`, `bias` (time axis) stay; `derived_update` computes `air_anomaly` after the atmosphere |
| `libs/yelmox_climate.f90` | step order: backend (snapshot terms) → uniform `tas`/`pr` terms → derived coefficients → ocean terms → taper → bias; contract checks at init |
| `yelmox_bipolar/obm_coupling.f90` | sets the `obm_*` derived coefficients from the boxes where `obm2ism` ran; `[ocn_obm_<sfx>]` removed |
| `libs/climate_rembo.f90` | reads the uniform `tas` term instead of `[anom]` |
| `input/yelmox_defaults.nml` | templates `coef`, `term`, `snapshot`, `snapshots`, `forcing`; the old groups removed |
| `input/forcing_db.nml` | `snapshot_*`, `coef_*`, `term_*` with descriptive names |
| `yelmox/*.nml`, `yelmox_bipolar/*.nml`, `scripts/*/params.sh` | converted; the `-p` arrays shrink to `forcing.terms` plus run control |
| `scripts/check_par_nml.py` | templates by prefix `coef_`, `term_`, `snapshot_` |
| `scripts/forcing-builder/`, `docs/forcing-builder/` | export the new groups; the builder rebuilt in the order of the equation (time, equation, pieces, after), the run's forcing groups as its state, the par file and `-p` as its export (done 2026-10-21) |
| docs | `yelmox.md`, `climate-snapshots.md`, experiment pages, CHANGELOG; "forcing time" wording |

## 8. Phases

All on dev, one branch `forcing-terms`, merged when the whole chain passes.

| # | Content | Check |
|---|---|---|
| 1 | this doc, the mock-up in `docs/` | |
| 2 | coefficient registry: the `derived` kind, `coefs_set` (done; `indices_*` renamed `coefs_*` on 2026-10-22, `tests/test_coefs.f90`) | unit: a derived coefficient set and read |
| 3 | `libs/climate_terms.f90`; terms in the snapshot backend, uniform and replace terms, derived `air_anomaly`, bridge sets the OBM coefficients; `anom`/`ocn_fraction`/`ocn_obm` removed; new keys, library, defaults, bases, bipolar and checker converted in the same commit (every commit runs) | `check_par_nml` clean, both drivers build |
| 4 | scripts (`-p` arrays, docs-figures, run_mg_forward, bipolar tests) | |
| 5 | docs-figures 18 cases + the bipolar case on Levante | bit-identical to dev 54d9ee6 (the same operations in the same order) |
| 6 | forcing builder on the new groups, equation first (the state the forcing groups, import from a base plus `-p`, export the difference); docs; CHANGELOG | done |

## 9. Decisions taken (2026-10-21)

Encoding (a), one group per term; (b) recorded in section 6. The uniform
anomaly is a term. Derived coefficients as a coefficient kind, with the factor in
the coefficient group. Fields always explicit. Precipitation by the fixed rule, nothing to
set. `snap0`/`snap1` for the two snapshots of a pattern (0 and 1 as the
coefficient values, tied to `norm_t0`/`norm_t1`); `f_ta/f_to/f_so` kept;
`tstep_method` kept; `forcing_time` wording, `clock` dropped. OBM absolute
mode as the `replace` pattern. Time-dependent bias through the varslice time
axis, not as a term. Normalisation windows kept (`norm_t0 = t` or `t_a t_b`).
Taper, bias, orbit, sea level, border ice stay outside the sum.

Open: the group names `[forcing]`, `[snapshots]`, prefixes `term_`,
`snapshot_` (as on the mock-up unless changed; `series_` became `coef_` on
2026-10-22). (`[snapshots] combine` had no use and was removed.)
