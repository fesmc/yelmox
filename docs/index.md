# YelmoX

**YelmoX** is the driver framework around the [Yelmo](https://github.com/fesmc/yelmo) ice-sheet model. It bundles Yelmo together with the libraries needed to run realistic ice-sheet simulations: climate and ocean forcing modules, the [FastIsostasy](https://github.com/palma-ice/FastIsostasy) bedrock module, the [REMBOv1](https://github.com/alex-robinson/rembo1) energy/moisture-balance model, and shared utilities from [fesm-utils](https://github.com/fesmc/fesm-utils).

The YelmoX code repository is here: [https://github.com/fesmc/yelmox](https://github.com/fesmc/yelmox).

For documentation of the Yelmo ice-sheet model itself, see [the Yelmo docs](https://fesmc.github.io/yelmo/).

## Super-quick start

A summary of commands to get YelmoX running is given below. The example uses the `dkrz_levante_ifx` config — substitute your own config file from `config/` for a different system. Where data paths (`ice_data`, `isostasy_data`) appear, replace `/path/to/...` with the location on your system.

```bash
# Pick a working directory where all the repos will live as siblings
mkdir -p ~/models && cd ~/models

# ---------------------------------------------------------------
# fesm-utils — provides LIS, FFTW, and shared utility modules.
# Needed by yelmox, yelmo, and FastIsostasy.
# ---------------------------------------------------------------
git clone git@github.com:fesmc/fesm-utils.git
cd fesm-utils
./install_dkrz.sh ifx                       # installs lis-serial, lis-omp, fftw-serial, fftw-omp
cd utils
python3 config.py config/dkrz_levante_ifx
make clean && make fesmutils-static
cd ..
export FESMUSRC=$PWD                        # absolute path to the fesm-utils root
cd ..

# ---------------------------------------------------------------
# coordinates — needed by rembo1
# ---------------------------------------------------------------
git clone git@github.com:fesmc/coordinates.git
cd coordinates
python3 config.py config/dkrz_levante_ifx
make clean && make coord-static
export COORDSRC=$PWD
cd ..

# ---------------------------------------------------------------
# yelmo
# ---------------------------------------------------------------
git clone git@github.com:fesmc/yelmo.git
cd yelmo
python3 config.py config/dkrz_levante_ifx
ln -s $FESMUSRC ./
cd ..

# ---------------------------------------------------------------
# FastIsostasy
# ---------------------------------------------------------------
git clone git@github.com:palma-ice/FastIsostasy.git
cd FastIsostasy
python3 config.py config/dkrz_levante_ifx
ln -s $FESMUSRC ./
cd ..

# ---------------------------------------------------------------
# rembo1
# ---------------------------------------------------------------
git clone git@github.com:alex-robinson/rembo1.git
cd rembo1
python3 config.py config/dkrz_levante_ifx
ln -s $COORDSRC libs/
cd ..

# ---------------------------------------------------------------
# yelmox itself
# ---------------------------------------------------------------
git clone git@github.com:fesmc/yelmox.git
cd yelmox
python3 config.py config/dkrz_levante_ifx

# Link the component repos and fesm-utils into the yelmox root
ln -s $FESMUSRC ./
ln -s ../yelmo ./
ln -s ../FastIsostasy ./
ln -s ../rembo1 ./

# Link to the input data repositories (replace with your local paths)
ln -s /path/to/ice_data
ln -s /path/to/isostasy_data

# Copy the runme config file to the main directory and edit for your system
cp .runme/runme_config .runme_config

# Compile the default yelmox program
make clean
make yelmox

# Run a test simulation of Antarctica for 1000 yrs
./runme -r -e yelmox -n par/yelmo_Antarctica.nml -o output/ant-test -p ctrl.time_end=1e3
```

That's it!

### Notes

- The `FESMUSRC` environment variable should be set to the **absolute** path of the `fesm-utils` root directory (the one containing `utils/`, `fftw-serial/`, `lis-serial/`, etc.), not the `utils/` subdir.
- `make yelmox` will recursively build the static libraries for `fesm-utils/utils`, `yelmo`, `FastIsostasy`, and link them. To build the REMBO-coupled variant instead, use `make yelmox_rembo` (this also pulls in `rembo1`).
- Other available programs include `yelmox_ismip6`, `yelmox_esm`, `yelmox_rtip`, `yelmox_nahosmip`, `yelmox_bipolar`, and `yelmox_glaciers` — see the `Makefile` for the full list.
- For HPC-specific notes (modules, environment variables, queue setup) see the [Yelmo HPC notes](https://fesmc.github.io/yelmo/hpc-notes.html).
