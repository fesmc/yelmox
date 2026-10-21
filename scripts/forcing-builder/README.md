# forcing-builder

The data and the command-line side of the forcing builder page
(`docs/forcing-builder.md`, published with the docs): a run's configuration is
its base par file plus runme `-p` overrides, resolved against the forcing
library (`input/forcing_db.nml`) and the schema defaults
(`input/yelmox_defaults.nml`; the forcing as terms: `docs/_forcing-terms.md`),
and the forcing it imposes is drawn on the run's timeline.

Files:

- `Namelist.jl`: reader of the namelist files, keeping the comments (the
  schema's help text).
- `TsGen.jl`: the series methods of fesm-utils `tsgen` as pure functions.
- `Resolve.jl`: base + library + defaults + overrides, with the source of every
  value (`.runme/info.json` `par_defaults`, `par_library`).
- `export.jl`: writes `docs/forcing-builder/forcing.json` (schema, library,
  the run-control and forcing groups of the bases, the climate database groups,
  the scalar series of `input/` downsampled, and the cases of the docs figures
  as presets, from `scripts/docs-figures/run.sh --list`).
- `explain.jl`: the same resolution on the command line (the forcing time,
  the series, the terms with their equations and snapshots, the transforms,
  the orbit, the sea level), with a CairoMakie figure.

## Export (after changing a par file, the library, the defaults, a series or a docs-figures case)

```bash
julia --project=scripts/forcing-builder -e 'using Pkg; Pkg.instantiate()'
```

```bash
julia --project=scripts/forcing-builder scripts/forcing-builder/export.jl
```

The JSON is committed with the docs, so the published page needs no Julia; run
the export before publishing the docs.

## Explain

```bash
julia --project=scripts/forcing-builder scripts/forcing-builder/explain.jl yelmox/yelmox_Greenland.nml \
      -p ctrl.tstep_method=rel ctrl.time_init=-120e3 ctrl.time_end=0 forcing.series==orbital,dTo,ramp forcing.terms==orbital,ocean,ramp
```

```bash
julia --project=scripts/forcing-builder scripts/forcing-builder/explain.jl output/grl-ramp
```

A par file with `-p` options, or a run folder (its staged par file, the figure
written next to it as `forcing.png`). `--fig OUT.png` names the figure, `--fig
none` skips it.
