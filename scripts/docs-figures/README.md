# docs-figures

The forcing figures of the experiment pages (`docs/figures/<case>.png`): each
test case is run climate-only (ice sheet off, present-day topography) over a
forcing time range, and the climate at chosen forcing times is mapped against a
reference time: annual and summer air temperature difference, precipitation
ratio, and the surface mass balance difference when the climate supplies it
(ESM), or the reference climate itself for present-day cases.

## Run (Levante)

```bash
scripts/docs-figures/run.sh
```

The cases (par file, timeline, overrides) are the table in `run.sh`; give case
names to run only those. Output: `output/docs-figures/<case>` (`OUT`). Slurm
options: `RUNOPTS` (default `-q compute -w 00:30:00 --omp 16`); `STAGE=1` stages
the runs only.

## Plot

```bash
julia --project=scripts/docs-figures -e 'using Pkg; Pkg.instantiate()'
```

```bash
scripts/docs-figures/plot.sh
```

The figures (forcing times, reference, title) are the table in `plot.sh`; give
figure names to draw only those. `plot.jl` draws one figure; see its header for
the options (`--ref`, `--at`, `--ref-run`, `--abs`, `--title`).
