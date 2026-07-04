#!/usr/bin/env bash
#
# Antarctica multigrid resolution test (yelmox).
#
# Optimized present-day climate (no warming), 25 kyr each. Six runs:
#   - three "core" runs: Yelmo and every module on the same grid
#       (ANT-8KM, ANT-16KM, ANT-32KM);
#   - three multigrid runs: marine_shelf on a grid FINER than Yelmo
#       (Yelmo 32KM + mshlf 8KM / 16KM, Yelmo 16KM + mshlf 8KM).
#
# All runs use equil_method=opt (the nml default): basal friction (cf_ref) and
# thermal forcing (tf_corr) are optimized toward the PD target over the default
# 0-15 kyr windows, then free-evolve to time_end=25 kyr as a relaxation tail.
# The pmpt grounding-zone bmb scaling uses H_t=100 (gz_Hg1=100). 2D output every
# 1000 yr, including the diagnosed initial state.
#
# Grid mapping per run: Yelmo, isostasy and smb share <yelmo_grid>; marine_shelf
# and the hi-res topography hub share <mshlf_grid> (the finest grid, >= Yelmo
# resolution). For a core run the two grids are equal. Climate (snapclim
# atmosphere + ocean) always runs on ANT-32KM -- the only grid with the ISMIP6
# ocean forcing -- and is remapped by the coupler to the consumer grids.
#
# Usage (from anywhere):
#     yelmox/run_mg_resolution.sh
#
# Submits to the queue by default. Set runopts='-r' to run locally instead.

cd "$(dirname "$0")/.." || exit 1        # repo root

runopts='-rs -q 12h -w 10:00:00'

EXE="yelmox"
NML="yelmox/yelmox_Antarctica.nml"
OUTROOT="output/mg_opt"

# run_case <yelmo_grid> <mshlf_grid>
run_case() {
    local ygrid="$1" mgrid="$2"
    local out="$OUTROOT/y${ygrid#ANT-}_m${mgrid#ANT-}"
    runme $runopts -e "$EXE" -n "$NML" -o "$out" \
        -p yelmo.grid_name="$ygrid" htopo.grid_name="$mgrid" \
           coupling.grid_mshlf="$mgrid" coupling.grid_isos="$ygrid" \
           coupling.grid_clim=ANT-32KM coupling.grid_smb="$ygrid" \
           ctrl.time_end=25e3 \
           tm_2D.dt=1000 ytopo.gz_Hg1=100
}

# --- core runs: Yelmo and every module on the same grid ---
run_case ANT-8KM  ANT-8KM
run_case ANT-16KM ANT-16KM
run_case ANT-32KM ANT-32KM

# --- marine_shelf on a finer grid than Yelmo ---
run_case ANT-32KM ANT-8KM
run_case ANT-32KM ANT-16KM
run_case ANT-16KM ANT-8KM
