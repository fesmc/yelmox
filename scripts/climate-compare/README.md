# climate-compare

Climate-only comparison of two configurations. Two par files (A and B) are run
with the ice sheet off and the climate updated every `dt` years over the forcing
range; the climate 2D files of both runs are then compared: domain means of
`t2m_ann`, `pr_ann` and `to_ann` per time (printed table), and maps of A, B and
B−A with the time series of the means (one figure per field).

It is the regression check for changes to the climate backends and modifiers
(`docs/_forcing-design.md` section 13). It is also where the open data questions
of the forcing refactor are answered: whether `clim_sigma` is still needed after
the global remap of the snapshots (section 9), and the I/O cost per update of
time-varying snapshots (GRL paleo; `timing.txt` of each run).

## Run (Levante)

```bash
scripts/climate-compare/run.sh --pb "coef_dTo.factor=0.5" \
    yelmox/yelmox_Greenland.nml yelmox/yelmox_Greenland.nml output/climcmp/grl
```

Options: `--t0`, `--t1` (forcing range, default -130e3 to 0), `--pd` (present
day: forcing time fixed at 0), `--dt` (default 1000 yr), `-p "k=v ..."`
(overrides for both), `--pa`/`--pb` (A or B only), `--stage` (submit scripts
only), `--local` (no slurm). Slurm options: `RUNOPTS` (default
`-q compute -w 00:30:00 --omp 16`). The jobs are named `climcmp-A`, `climcmp-B`.

Overrides set on both runs: `comps.with_ice_sheet=False`, `comps.dt_clim=dt`,
the timeline (`tstep_method = "rel"`, `time_init`, `time_end`, `dtt = dt`; in the
group named by `ctrl.run_step`), the output cadence (`tm_2D`, `tm_2Dsm`, `tm_1D`
every `dt`) and no restart bundles but the final one.

## Compare

```bash
julia --project=scripts/climate-compare -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=scripts/climate-compare scripts/climate-compare/compare.jl output/climcmp/grl
```

Options: `--file` (default: `snap.nc`, `esm.nc` or `rembo.nc`, whichever is in
`A/`), `--depth D` (`to_ann` at the level nearest `D` m; default the depth mean),
`--ext pdf`. Figures: `OUTDIR/climcmp_<field>.png`. The "ice" means use
`H_ice > 0` of the first record of `A/yelmo.nc` (the reference topography, the
ice sheet being off) when it is on the climate grid; otherwise all cells. A field
missing in A or B is skipped.
