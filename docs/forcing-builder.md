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
terms](yelmox.md#forcing)): a reference snapshot plus a sum of terms, each a
coefficient times a pattern on explicit fields. Left, top to bottom: the
forcing time τ; the equation, first the reference, then the terms, then the
coefficients they name (a constant, a record, a ramp, a controller, or derived
from the model), each in its own card with where its data comes from (the
record of a coefficient; the clim_db fields of each snapshot a term names) and
what is missing; and what comes after the composition. Right: the equation
chain evaluated at a chosen time, the coefficients over the run, and the run
written out as a par file and as the `-p` list with the `runme` command against
a base par file. Start from a case of the docs figures or from a base par file
as it is.

<details class="fb-about"><summary>How it is built, and its limits</summary>

The page holds the forcing groups of the run (`[forcing]`, `[series_*]`,
`[term_*]`, `[snapshot_*]`, `[snapshots]`, the transforms, the time group).
The page says "coefficient" for the scalar α(τ) of a term; in the par file it
stays a group `[series_<name>]`, listed in `[forcing] series` and named by
the key `series` of a term.
A group is taken from the base par file (`yelmox/yelmox_<domain>.nml`), else
from the forcing library (`input/forcing_db.nml`), else created new from its
template in `input/yelmox_defaults.nml` when the equation names it, and edited
in place; writing the run out compares each group with its source, so the
`-p` list is what differs and a new group is given as a block to add to the
par file or the library. A preset or a base plus `-p` is read the other way
into the groups. The equation is built terms first: a term names its
coefficient and its snapshots from the ones defined, or a new name, and the
right column shows the whole composition at once; a coefficient a term names
but the equation does not define yet is marked, with a button that defines it
(a constant 1 to start from). A symbol is complete where it is introduced: a
record coefficient shows its file, normalisation and the record over its full
extent; a snapshot term, and the reference, show each snapshot with the
`input/clim_db.nml` groups of every field, marked where the equation needs it,
and a table of what each group reads (variable, file, τ range; its entry on a
click). A snapshot named by two terms is the same group in both cards.
The files are read as exported to `forcing-builder/forcing.json`
by `scripts/forcing-builder/export.jl`; the page has no server. The rare keys
of a group (noise, windows, restarts, smoothing) are kept as the par file sets
them and listed, not edited. The coefficient curves follow the fesm-utils `tsgen`
formulas; `scripts/forcing-builder/explain.jl` prints the same resolution and
draws the same figures on the command line, for a par file with `-p` options
or an existing run folder.

Limits: scalar coefficients only (the 2D snapshot fields are named and
checked, not drawn; a derived coefficient, set by the model, is shown by its
definition); the
ESM backend is one block (it composes its own climate, the terms add to it);
single-domain `yelmox` groups (no domain suffix); the `-p` list is a fixed
override per key, not an ensemble.

</details>

<div id="fb"><p>Loading…</p></div>

:::

<script type="module" src="forcing-builder/app.js"></script>
