#!/usr/bin/env bash
#
# ISMIP7 workflow for Greenland (yelmox, climate = esm).
#
#   Step 1  spinup     15-kyr present-day OPTIMIZED ice-sheet spin-up
#                      (sim.opt=True) -> writes a restart bundle
#   Step 2  scenarios  ssp126 / ssp370 / ssp585, each branched off that bundle
#
# The ice sheet + isostasy are ACTIVE (comps.with_ice_sheet/with_isostasy=True in
# yelmox_esm_Greenland.nml). The spin-up optimizes basal friction + thermal forcing to
# present day (&opt cf/tf_time_end=15e3), and the scenarios evolve the ice sheet under
# ISMIP7 climate/ocean forcing.
#
# Run the steps in order on the cluster; let the spin-up finish before launching
# the scenarios (they read its restart bundle):
#
#   scripts/ismip7/run_ismip7_greenland.sh spinup
#   scripts/ismip7/run_ismip7_greenland.sh scenarios
#
# Stage only -- create the run dirs + SLURM submit script but do NOT submit
# (a dry run to inspect everything first): set STAGE=1
#
#   STAGE=1 scripts/ismip7/run_ismip7_greenland.sh spinup
#   STAGE=1 scripts/ismip7/run_ismip7_greenland.sh scenarios
#
set -euo pipefail
cd "$(dirname "$0")/../.." || exit 1               # repo root

# ---- configuration ---------------------------------------------------------
EXE="yelmox"                                       # -> libyelmox/bin/yelmox.x (climate = esm)
NML="yelmox/yelmox_esm_Greenland.nml"
OUTROOT="output/ismip7_grl"
GCM="CESM2-WACCM"
GRID="GRL-8KM"
SCENARIOS=(ssp126 ssp370 ssp585)

SPINUP_YEARS=15000                                 # ice-sheet opt spin-up (matches &opt cf/tf_time_end=15e3)
PROJ_END=2300                                      # scenario end year (CE)

# runme submit options. STAGE=1 writes the submit script without submitting.
if [ "${STAGE:-0}" = 1 ]; then SUBMIT="-s"; else SUBMIT="-rs"; fi
# Walltimes differ hugely: the 15-kyr opt spin-up is the long job; the ~285-yr
# projections are short. Spin-up estimate from run_mg_resolution.sh (25 kyr GRL ~10h);
# verify against your actual throughput + queue max walltime (bump to a longer queue
# if 8 h is not enough).
HPCOPT_SPINUP="-q compute -w 08:00:00 --omp 8"
HPCOPT_SCEN="-q compute -w 02:00:00 --omp 8"

SPINUP_OUT="$OUTROOT/spinup"
# The spin-up's final restart bundle. yelmox names it restart-<time/1e3 %.3f>-kyr.
# Absolute path: the executable runs from inside the scenario's run dir, so a
# repo-root-relative path would not resolve.
BUNDLE="$(pwd)/$SPINUP_OUT/restart-$(awk "BEGIN{printf \"%.3f\", $SPINUP_YEARS/1000}")-kyr"

# ---- steps -----------------------------------------------------------------
case "${1:-}" in
  spinup)
    runme $SUBMIT $HPCOPT_SPINUP -e "$EXE" -n "$NML" -o "$SPINUP_OUT" \
      -p ctrl.run_step=spinup sim.relax=True sim.opt=True \
         domain.grid_hub="$GRID" \
         spinup.time_init=0 spinup.time_end="$SPINUP_YEARS"
    ;;
  scenarios)
    for exp in "${SCENARIOS[@]}"; do
      runme $SUBMIT $HPCOPT_SCEN -e "$EXE" -n "$NML" -o "$OUTROOT/$exp" \
        -p ctrl.run_step=transient esm.experiment="$exp" esm.esm_name="$GCM" \
           esm.use_esm=True esm.use_hist=True esm.use_proj=True \
           domain.grid_hub="$GRID" \
           ctrl.restart="$BUNDLE" \
           transient.time_init=2015 transient.time_end="$PROJ_END"
    done
    ;;
  *)
    echo "usage: $0 {spinup|scenarios}   (prefix STAGE=1 to stage without submitting)" >&2
    exit 1
    ;;
esac
