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

Start from a preset or a base par file and change the run on the left; the
right column writes the forcing as equations at a chosen forcing time τ,
draws the series over the run, and writes the run out three ways: the par
file, the equation as text, and the `-p` list with the `runme` command.

<details class="fb-about"><summary>How it is built, and its limits</summary>

A run's configuration is its base par file (`yelmox/yelmox_<domain>.nml`) plus
`-p` overrides, resolved against the forcing library (`input/forcing_db.nml`)
and the schema defaults (`input/yelmox_defaults.nml`), see
[yelmox](yelmox.md#forcing). The page reads those files as exported to
`forcing-builder/forcing.json` by `scripts/forcing-builder/export.jl`; it has
no server. Every widget is one namelist key, labelled with where its value
comes from (base, library, default, override). The presets are the cases of
the docs figures (`scripts/docs-figures/run.sh`). The series curves follow the
fesm-utils `tsgen` formulas; `scripts/forcing-builder/explain.jl` prints the
same resolution and draws the same figures on the command line, for a par
file with `-p` options or an existing run folder.

Limits: scalar series only (the 2D snapshot fields are listed, not drawn; a
derived series, set by the model, is shown by its definition); single-domain
`yelmox` groups (no domain suffix); the `-p` list is a fixed override per key,
not an ensemble.

</details>

<div id="fb"><p>Loading…</p></div>

:::

<script type="module" src="forcing-builder/app.js"></script>
