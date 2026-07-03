#!/usr/bin/env bash
#
# Antarctica multigrid resolution test (yelmox_mg).
#
# Steady-state, present-day climate (no warming), 1000 yr each. Six runs:
#   - three "core" runs: Yelmo and every module on the same grid
#       (ANT-8KM, ANT-16KM, ANT-32KM);
#   - three multigrid runs: marine_shelf on a grid FINER than Yelmo
#       (Yelmo 32KM + mshlf 8KM / 16KM, Yelmo 16KM + mshlf 8KM).
#
# All runs use equil_method=none (free evolution) and the pmpt grounding-zone
# bmb scaling with H_t=100 (gz_Hg1=100). 2D output every 100 yr, including the
# diagnosed initial state.
#
# Grid mapping per run: Yelmo, isostasy, climate and smb share <yelmo_grid>;
# marine_shelf and the hi-res topography hub share <mshlf_grid> (the finest
# grid, >= Yelmo resolution). For a core run the two grids are equal.
#
# Usage (from anywhere):
#     yelmox_mg/run_mg_resolution.sh
#
# Submits to the queue by default. Set runopts='-r' to run locally instead.

cd "$(dirname "$0")/.." || exit 1        # repo root

runopts='-rs -q 12h -w 06:00:00'

EXE="mg"
NML="yelmox_mg/yelmox_mg_Antarctica.nml"
OUTROOT="output/mg"

# run_case <yelmo_grid> <mshlf_grid>
run_case() {
    local ygrid="$1" mgrid="$2"
    local out="$OUTROOT/y${ygrid#ANT-}_m${mgrid#ANT-}"
    runme $runopts -e "$EXE" -n "$NML" -o "$out" \
        -p yelmo.grid_name="$ygrid" htopo.grid_name="$mgrid" \
           coupling.grid_mshlf="$mgrid" coupling.grid_isos="$ygrid" \
           coupling.grid_clim="$ygrid" coupling.grid_smb="$ygrid" \
           coupling.equil_method=none ctrl.time_end=1000 \
           tm_2D.dt=100 ytopo.gz_Hg1=100
}

# --- core runs: Yelmo and every module on the same grid ---
run_case ANT-8KM  ANT-8KM
run_case ANT-16KM ANT-16KM
run_case ANT-32KM ANT-32KM

# --- marine_shelf on a finer grid than Yelmo ---
run_case ANT-32KM ANT-8KM
run_case ANT-32KM ANT-16KM
run_case ANT-16KM ANT-8KM
