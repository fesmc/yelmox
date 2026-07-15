#!/usr/bin/env bash
#
# Run LIS smb_simple simulations at different times (LGM, MIS-3) using the
# Batchelor2019 target ice mask, and optionally build climate/ice-sheet
# parameter ensembles.
#
# Requires v2.3+ (var_io isochrone table fix). Jobs are submitted via SLURM
# (`runme -rs`) rather than run on the login node, which gets killed for
# longer runs.
#
# Usage:
#   scripts/run_lis_smbsimple.sh              # submit the two short single tests
#   RUN_ENSEMBLES=1 scripts/run_lis_smbsimple.sh   # also submit the ensembles
#
set -euo pipefail

# ---- configuration -----------------------------------------------------------
NML=yelmox/yelmox_LIS.nml                 # LIS namelist (smb_method="smb_simple" default)
EXE=yelmox                                # runme exe alias -> libyelmox/bin/yelmox.x
RUNOPT="-rs -q 12h -w 10:00:00"           # SLURM submit; tune queue/walltime
OUTROOT="tmp/lischeck/$(date +%Y-%m-%d)"  # base output dir
TIME_END=30e3                             # short test length; raise (e.g. 30e3) for production
RUN_ENSEMBLES="${RUN_ENSEMBLES:-0}"       # set to 1 to also submit the parameter ensembles

# Batchelor2019 mask index -> label (variable `label` in
# ice_data/Laurentide/LIS-32KM/LIS-32KM_Batchelor2019_ice_masks.nc):
#   14=45ka  15=40ka  16=35ka (MIS-3)  17=30ka  18=LGM
# ctrl.tstep_const sets the assumed climate time (orbit/insolation), in yr BP.

# ---- helper: submit one single run or ensemble for a given period ------------
run_period () {   # args: name  mask_idx  tstep_const  [extra -p key=val ...]
    local name=$1 mask_idx=$2 tstep=$3
    shift 3
    echo ">> $name  (mask_idx=$mask_idx, tstep_const=$tstep, time_end=$TIME_END)"
    runme $RUNOPT -e "$EXE" -n "$NML" -o "$OUTROOT/$name" \
        -p ctrl.time_end="$TIME_END" ctrl.tstep_const="$tstep" smb_simple.mask_idx="$mask_idx" "$@"
}

# ---- 1) single short test runs ----------------------------------------------
run_period lgm-test  18 -21e3
run_period mis3-test 16 -35e3

# ---- 2) parameter ensembles --------------------------------------------------
# Each comma list (a=1,2,3) adds an ensemble dimension; multiple lists form a
# full (Cartesian) grid of members under the -o directory.
#
# Ice sheet: basal friction (ytill) + effective pressure (yhyd).
#   NOTE: legacy `yneff.delta` is now `yhyd.till_delta`.
ICE_PARAMS="ytill.cf_ref=0.05,0.1,0.2 ytill.scale_zb=2 ytill.z0=-100 ytill.z1=200 yhyd.till_delta=0.1"

# Climate: smb_simple snowline intercept a1 is the main knob (m);
#   see also smb_simple.beta0 (ablation lapse) and smb_simple.c0 (accumulation).
CLIM_PARAMS="smb_simple.a1=3700,3900,4100"

if [ "$RUN_ENSEMBLES" = "1" ]; then
    run_period lgm-ens  18 -21e3 $ICE_PARAMS $CLIM_PARAMS   # 3 cf_ref x 3 a1 = 9 members
    run_period mis3-ens 16 -35e3 $ICE_PARAMS $CLIM_PARAMS
fi

echo "done."
