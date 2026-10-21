---
title: "Forcing builder"
page-layout: custom
toc: false
format:
  html:
    include-in-header:
      text: |
        <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
        <link rel="stylesheet" href="forcing-builder/app.css">
---

::: {.fb-page}

The forcing of a run, built in the order of its equation (the [forcing as
terms](yelmox.md#forcing)). Left, top to bottom: the forcing time τ; the
equation, a reference snapshot plus terms, each a scalar series times a
pattern on explicit fields; the pieces the equation names, the record of each
series and the fields of each snapshot, with what is missing; and what comes
after the composition. Right: the equation chain evaluated at a chosen time,
the series over the run, and the run written out as a par file and as the
`-p` list with the `runme` command against a base par file. Start from a case
of the docs figures or from a base par file as it is.

<details class="fb-about"><summary>How it is built, and its limits</summary>

The page holds the forcing groups of the run (`[forcing]`, `[series_*]`,
`[term_*]`, `[snapshot_*]`, `[snapshots]`, the transforms, the time group).
A group is taken from the base par file (`yelmox/yelmox_<domain>.nml`), else
from the forcing library (`input/forcing_db.nml`), else created new from its
template in `input/yelmox_defaults.nml` when the equation names it, and edited
in place; writing the run out compares each group with its source, so the
`-p` list is what differs and a new group is given as a block to add to the
par file or the library. A preset or a base plus `-p` is read the other way
into the groups. The files are read as exported to `forcing-builder/forcing.json`
by `scripts/forcing-builder/export.jl`; the page has no server. The rare keys
of a group (noise, windows, restarts, smoothing) are kept as the par file sets
them and listed, not edited. The series curves follow the fesm-utils `tsgen`
formulas; `scripts/forcing-builder/explain.jl` prints the same resolution and
draws the same figures on the command line, for a par file with `-p` options
or an existing run folder.

Limits: scalar series only (the 2D snapshot fields are named and checked, not
drawn; a derived series, set by the model, is shown by its definition); the
ESM backend is one block (it composes its own climate, the terms add to it);
single-domain `yelmox` groups (no domain suffix); the `-p` list is a fixed
override per key, not an ensemble.

</details>

<div id="fb"><p>Loading…</p></div>

:::

<script type="module" src="forcing-builder/app.js"></script>
