# YelmoX

[![Docs](https://img.shields.io/badge/docs-quarto-75aadb?logo=quarto&logoColor=white)](https://fesmc.github.io/yelmox/)

**YelmoX** is the driver framework around the [Yelmo](https://github.com/fesmc/yelmo) ice-sheet model ([docs](https://fesmc.github.io/yelmo/)). It bundles Yelmo together with the libraries needed to run realistic ice-sheet simulations: climate and ocean forcing modules, the [FastIsostasy](https://github.com/palma-ice/FastIsostasy) bedrock module, the [REMBOv1](https://github.com/alex-robinson/rembo1) energy/moisture-balance model, and shared utilities from [fesm-utils](https://github.com/fesmc/fesm-utils).

See the [YelmoX documentation](https://fesmc.github.io/yelmox/) for installation and usage.

> **Compatibility (dev): fesm-utils API shift.** YelmoX now requires **fesm-utils dev at
> `3f415cc` (2026-06-26) or later**. fesm-utils folded its standalone `coordinates` library
> into `utils/src/coords/` and moved `mv`/`TOL` from `precision` to a new `constants` module
> and `nc_read_interp` from `mapping_scrip` to `ncio_interp` (`mps=(map_scrip_class)` →
> `map=(map_class)`). `ismip6` and `marine_shelf` were updated accordingly, and `rembo1` now
> relies on fesm-utils for its shared utility modules. All bundled libraries (yelmo,
> FastIsostasy, FastHydrology, rembo1) must be rebuilt against the same fesm-utils.
