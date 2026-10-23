# YelmoX

[![Docs](https://img.shields.io/badge/docs-quarto-75aadb?logo=quarto&logoColor=white)](https://fesmc.github.io/yelmox/)

**YelmoX** is the driver framework around the [Yelmo](https://github.com/fesmc/yelmo) ice-sheet model ([docs](https://fesmc.github.io/yelmo/)). It bundles Yelmo together with the libraries needed to run realistic ice-sheet simulations: climate and ocean forcing modules, the [chion](https://github.com/fesmc/chion) surface mass balance model, the [FastIsostasy](https://github.com/palma-ice/FastIsostasy) bedrock module, the [REMBOv1](https://github.com/alex-robinson/rembo1) energy/moisture-balance model, and shared utilities from [fesm-utils](https://github.com/fesmc/fesm-utils).

See the [YelmoX documentation](https://fesmc.github.io/yelmox/) for installation and usage.

> **Dependencies (v2.4).** YelmoX builds against
> [fesm-utils](https://github.com/fesmc/fesm-utils),
> [yelmo](https://github.com/fesmc/yelmo) (with
> [FastHydrology](https://github.com/fesmc/FastHydrology), cloned inside it),
> [FastIsostasy](https://github.com/palma-ice/FastIsostasy),
> [rembo1](https://github.com/alex-robinson/rembo1),
> [chion](https://github.com/fesmc/chion) (surface mass balance) and
> [insol](https://github.com/fesmc/insol) (insolation). The versions are pinned in
> `.configme/manifest.toml`: on main/stable the `stable` branch of each
> dependency, at a tagged release the dependency tags. Install either the latest
> stable versions (`configme install yelmox`) or a fixed tagged release
> (`configme install yelmox:v2.4`). All libraries must be built against the same
> fesm-utils.
