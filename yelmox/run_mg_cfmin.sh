#!/usr/bin/env bash
#
# Antarctica multigrid resolution test -- Yelmo-32KM reruns with a lower basal
# friction floor (yelmox).
#
# Reruns the three Yelmo-32KM spin-ups from run_mg_resolution.sh (the ones that
# completed) with the optimization's basal-friction lower bound reduced from the
# default 0.01 to 1e-4. That bound is set by [ytill] cf_min, because [opt]
# use_yelmo_cf_min=.true. routes the opt lower bound to Yelmo's till cf_min (the
# [opt] cf_min / opt_cf_min are bypassed while that switch is on). A floor of
# 0.01 was too high -- it prevented the optimizer from making ice-stream beds
# slippery enough, leaving the spun-up ice too thick. cf_min is passed as a
# parameter override, so the shared namelist is left untouched.
#
# Everything else matches run_mg_resolution.sh (opt spin-up, 25 kyr, same grids,
# 1000 yr 2D output). Output goes to a NEW folder so the original mg_opt runs are
# preserved for comparison.
#
# Usage (from anywhere):
#     yelmox/run_mg_cfmin.sh
#
# Submits to the queue by default. Set runopts='-r' to run locally instead.

cd "$(dirname "$0")/.." || exit 1        # repo root

runopts='-rs -q compute -w 08:00:00'

EXE="yelmox"
NML="yelmox/yelmox_Antarctica.nml"
OUTROOT="output/mg_opt_cfmin1e-4"        # new folder; original output/mg_opt kept
CFMIN=1e-4                               # basal-friction floor (was 0.01)

# run_case <yelmo_grid> <mshlf_grid>
run_case() {
    local ygrid="$1" mgrid="$2"
    local out="$OUTROOT/y${ygrid#ANT-}_m${mgrid#ANT-}"
    runme $runopts -e "$EXE" -n "$NML" -o "$out" \
        -p yelmo.grid_name="$ygrid" htopo.grid_name="$mgrid" \
           coupling.grid_mshlf="$mgrid" coupling.grid_isos="$ygrid" \
           coupling.grid_clim=ANT-32KM coupling.grid_smb="$ygrid" \
           ctrl.time_end=25e3 \
           tm_2D.dt=1000 ytopo.gz_Hg1=100 \
           ytill.cf_min="$CFMIN"
}

# --- the three Yelmo-32KM runs (core + two multigrid) ---
run_case ANT-32KM ANT-32KM
run_case ANT-32KM ANT-16KM
run_case ANT-32KM ANT-8KM
