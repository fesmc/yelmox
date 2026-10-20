---
title: "Forcing builder"
format:
  html:
    include-in-header:
      text: |
        <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
        <link rel="stylesheet" href="forcing-builder/app.css">
---

Choose a base par file, change the run control and the forcing, and read off
what the run will impose: the forcing clock, each index on the run's timeline,
the climate modifiers' anomalies, the sea level, with the equations behind them
and the values they take. The result is the `runme` command, with the `-p`
overrides that reproduce the choices, to copy into a script.

How it is built: a run's configuration is its base par file
(`yelmox/yelmox_<domain>.nml`) plus `-p` overrides, resolved against the forcing
library (`input/forcing_db.nml`) and the schema defaults
(`input/yelmox_defaults.nml`), see [yelmox](yelmox.md#forcing-library). The page
reads those files as exported to `forcing-builder/forcing.json` by
`scripts/forcing-builder/export.jl`; it has no server. Every widget is one
namelist key, labelled with where its value comes from (base, library,
default, override). The index curves follow the fesm-utils `tsgen` formulas;
`scripts/forcing-builder/explain.jl` prints the same resolution and draws the
same figures on the command line, for a par file with `-p` options or an
existing run folder.

Limits: scalar series only (the 2D snapshot fields are listed, not drawn);
single-domain `yelmox` groups (no domain suffix); the `-p` list is a fixed
override per key, not an ensemble.

<div id="fb"><p>Loading…</p></div>
<script type="module" src="forcing-builder/app.js"></script>
