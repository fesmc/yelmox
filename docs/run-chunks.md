---
title: "Long runs in chunks"
---

A long simulation can run as a chain of short jobs ("chunks") instead of one
long job. On a busy cluster short jobs start much sooner (the scheduler fits
them into gaps between larger jobs), and a chain of 30-min jobs often finishes
long before one 8-hour job starts. The restarts are exact, so the chain gives
the same result as one run.

- **Script:** `scripts/run_chunks.sh` (Slurm; stages every chunk with runme and
  submits it with a dependency on the previous one).

## How a chunk restarts

The start time of a run comes from its timeline group (`[ctrl]`, or the group
named by `[ctrl] run_step`, e.g. `[transient]`), not from the restart bundle.
Each chunk therefore gets

- `time_init` = the `time_end` of the previous chunk, and its own `time_end`;
- `ctrl.restart` = the bundle the previous chunk wrote at its `time_end`. A run
  always writes a bundle at `time_end`, named
  `restart-<time_end/1000, 3 decimals>-kyr` (e.g. `restart-15.000-kyr`,
  `restart--125.000-kyr`), so all paths are known when the chain is
  submitted.

Chunk k runs in its own folder `OUTDIR/cNN`; the output of the run is split
over the chunk folders. Each job waits for the previous one (`sbatch
--dependency=afterok`): if a chunk fails or runs out of time, the rest of the
chain is not started (Slurm keeps it pending with reason
`DependencyNeverSatisfied`; cancel it with `scancel`).

Choose the chunk length so a chunk finishes well within the walltime
(`RUNOPTS`, default 30 min with 16 OpenMP threads); the chunk boundaries should
be multiples of the climate, isostasy and output intervals (`dt_clim`,
`dt_isos`, `tm_*`), which are counted from the start of each run.

## What changes after the first chunk

Some settings are timed from the start of a run (`time_elapsed`), not from the
forcing clock, and must be set for the later chunks (`--p-rest`):

| Setting | Later chunks |
|---|---|
| Relaxation (`[sim] relax`, schedule `[relax] time1`, `time2`) | `sim.relax=False` once the relaxation has ended (let the first chunk cover it) |
| Optimization (`[opt] cf_time_*`, `tf_time_*`) | same: within the first chunk, or `sim.opt=False` after |
| elsa isochrones (`[ytrc] use_elsa`) | `ytrc.elsa_restart=True` (restore the layers) |
| Cold-start settings (`[sim] init_method`, ...) | none: a restart skips them |

`[yelmo] restart_z_bed = True` (all configs) takes the bed from the bundle; with
`False`, `restart_relax` would blend the bed towards the input topography after
each restart.

## Usage

```bash
scripts/run_chunks.sh [options] PAR OUTDIR T0 T1 [T2 ...]
```

`T0 T1 ...` are the chunk boundaries in the time of the timeline group (relative
time, calendar years or elapsed time, as the par file runs). Options:

- `-r BUNDLE`: restart bundle of the first chunk;
- `-d JOBID`: the first chunk waits for this job;
- `-p "k=v ..."`: overrides for all chunks; `--p-first`, `--p-rest`: for the first
  or the later chunks only;
- `-g GROUP`: the timeline group (default: `[ctrl] run_step`);
- `-e EXE`, `-n NAME`: runme executable alias, job-name prefix.

The script prints the job of each chunk and, on the last line, the job id and the
final bundle of the chain, so a next stage can be chained to it with `-d` and
`-r`.

## Example: Greenland paleo

The two [Greenland paleo](experiment-grl-paleo.md) setups, each a spin-up and a
transient, in 5-kyr chunks (the A spin-up in 6, 5 and 4 kyr, so that its first
chunk covers the relaxation until 6 kyr):

```bash
read jA rA < <(scripts/run_chunks.sh --p-rest "sim.relax=False" yelmox/yelmox_GreenlandPaleo_spinup.nml output/grl-paleo/spinup 0 6000 11000 15000)
```

```bash
scripts/run_chunks.sh -d "$jA" -r "$rA" --p-rest "ytrc.elsa_restart=True" yelmox/yelmox_GreenlandPaleo_lgp.nml output/grl-paleo/lgp $(seq -130000 5000 0)
```

```bash
read jB rB < <(scripts/run_chunks.sh yelmox/yelmox_GreenlandPaleo_lgm.nml output/grl-paleo/lgm 0 5000 10000 15000 20000)
```

```bash
scripts/run_chunks.sh -d "$jB" -r "$rB" --p-rest "ytrc.elsa_restart=True" yelmox/yelmox_GreenlandPaleo_deglac.nml output/grl-paleo/deglac -21000 -16000 -11000 -6000 -1000 0
```

The first command chain is 3 + 26 jobs, the second 4 + 5; all 38 are submitted
at once.
