# YelmoX

**YelmoX** is the driver framework around the [Yelmo](https://github.com/fesmc/yelmo) ice-sheet model. It bundles Yelmo together with the libraries needed to run realistic ice-sheet simulations: climate and ocean forcing modules, the [FastIsostasy](https://github.com/palma-ice/FastIsostasy) bedrock module, the [REMBOv1](https://github.com/alex-robinson/rembo1) energy/moisture-balance model, and shared utilities from [fesm-utils](https://github.com/fesmc/fesm-utils).

For documentation of the Yelmo ice-sheet model itself, see [the Yelmo docs](https://fesmc.github.io/yelmo/).

## Quick start

With [`configme`](https://github.com/fesmc/configme) installed (see [Install configme](#install-configme-one-time) below), build the whole YelmoX stack with a single command from the directory where you want the checkout to live:

```bash
configme install yelmox
```

That's it, you should now be ready to compile and run any yelmox program flavor:

```bash
make clean
make yelmox
```

## What `configme install yelmox` does

`configme install yelmox` clones YelmoX and its component repositories (`fesm-utils`, `coordinates`, `yelmo`, `FastIsostasy`, `rembo1`), configures each for your machine and compiler, links them into the YelmoX directory, and builds `fesm-utils` (LIS + FFTW + utils, which can take 10-30 min). If `configme` can detect your machine from the hostname it does so, otherwise it prompts you.

The only system dependency you must install yourself is **netCDF** (see the [Yelmo dependencies notes](https://fesmc.github.io/yelmo/getting-started.html#dependencies) for installation tips). Everything else — LIS, FFTW, the `fesm-utils` libraries, the component repositories, their links, and `runme` — is managed by `configme`.

Common options:

```bash
configme install yelmox -m dkrz_levante -c ifx   # pick the machine + compiler explicitly
configme install yelmox -d https                 # clone over HTTPS (no GitHub SSH key needed)
configme install yelmox --dir ~/models/yelmox    # put the checkout here instead of ./yelmox
configme install yelmox --overwrite              # re-clone over an existing checkout
configme install yelmox --build-deps             # rebuild dependency packages without prompting
```

Run `configme list` for the supported machines and compilers, and `configme --help` for the full command surface. The exact clone/configure/link/build commands `configme install yelmox` runs for you are recorded in a `.install.sh` script in the checkout, and are also shown for context on the [configme install details](configme-install-details.md) page — these are shown for reference only; `configme install` is the recommended path and you do not need to run them by hand.

## Manual clone of YelmoX

`configme install yelmox` already clones the repository for you, so most users do not need this step. If you want to fetch the YelmoX source tree by hand (for example, to inspect it before installing), the repository is here: [https://github.com/fesmc/yelmox](https://github.com/fesmc/yelmox).

```bash
git clone git@github.com:fesmc/yelmox.git
cd yelmox
```

## Install configme (one time)

YelmoX is configured and built with [`configme`](https://github.com/fesmc/configme), a small Python tool that detects your netCDF installation, configures every package in the stack for your machine and compiler, and clones/links/builds the whole thing with one command. It is installed once, globally, and provides the `configme` command on your `PATH`:

```bash
pip install git+https://github.com/fesmc/configme
```

To upgrade it later, use its built-in update command:

```bash
configme update
```

If the `configme` command is not found after installation, your Python user bin directory is probably not on your `PATH`; add it in your `~/.bashrc` / `~/.zshrc`:

```bash
export PATH="${PATH}:${HOME}/.local/bin"
```
