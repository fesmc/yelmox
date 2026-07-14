#!/usr/bin/env bash
#
# ISMIP7 forcing-only workflow for Greenland (refactored yelmox_esm).
#
#   Step 1  spinup     reference-climate spin-up  -> writes a restart bundle
#   Step 2  scenarios  ssp126 / ssp370 / ssp585, each branched off that bundle
#
# All runs are forcing-only (coupling.with_ice_sheet=False, with_isostasy=False):
# they produce ISMIP7 climate + ocean forcing, no ice dynamics.
#
# Run the steps in order on the cluster; let the spin-up finish before launching
# the scenarios (they read its restart bundle):
#
#   yelmox_esm/run_ismip7_greenland.sh spinup
#   yelmox_esm/run_ismip7_greenland.sh scenarios
#
# Stage only -- create the run dirs + SLURM submit script but do NOT submit
# (a dry run to inspect everything first): set STAGE=1
#
#   STAGE=1 yelmox_esm/run_ismip7_greenland.sh spinup
#   STAGE=1 yelmox_esm/run_ismip7_greenland.sh scenarios
#
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1                  # repo root

# ---- configuration ---------------------------------------------------------
EXE="esm"                                          # -> libyelmox/bin/yelmox_esm.x
NML="yelmox_esm/yelmox_esm_Greenland.nml"
OUTROOT="output/ismip7_grl"
GCM="CESM2-WACCM"
GRID="GRL-8KM"
SCENARIOS=(ssp126 ssp370 ssp585)

SPINUP_YEARS=10                                    # forcing-only: short suffices
PROJ_END=2300                                      # scenario end year (CE)

# runme submit options. STAGE=1 writes the submit script without submitting.
if [ "${STAGE:-0}" = 1 ]; then SUBMIT="-s"; else SUBMIT="-rs"; fi
HPCOPT="-q compute -w 02:00:00 --omp 8"

SPINUP_OUT="$OUTROOT/spinup"
# The spin-up's final restart bundle. yelmox names it restart-<time/1e3 %.3f>-kyr.
# Absolute path: the executable runs from inside the scenario's run dir, so a
# repo-root-relative path would not resolve.
BUNDLE="$(pwd)/$SPINUP_OUT/restart-$(awk "BEGIN{printf \"%.3f\", $SPINUP_YEARS/1000}")-kyr"

# ---- steps -----------------------------------------------------------------
case "${1:-}" in
  spinup)
    runme $SUBMIT $HPCOPT -e "$EXE" -n "$NML" -o "$SPINUP_OUT" \
      -p ctrl.run_step=spinup coupling.equil_method="opt" \
         yelmo.grid_name="$GRID" htopo.grid_name="$GRID" \
         spinup.time_init=0 spinup.time_end="$SPINUP_YEARS"
    ;;
  scenarios)
    for exp in "${SCENARIOS[@]}"; do
      runme $SUBMIT $HPCOPT -e "$EXE" -n "$NML" -o "$OUTROOT/$exp" \
        -p ctrl.run_step=transient esm.experiment="$exp" esm.esm_name="$GCM" \
           esm.use_esm=True esm.use_hist=True esm.use_proj=True \
           yelmo.grid_name="$GRID" htopo.grid_name="$GRID" \
           coupling.restart="$BUNDLE" \
           transient.time_init=2015 transient.time_end="$PROJ_END"
    done
    ;;
  *)
    echo "usage: $0 {spinup|scenarios}   (prefix STAGE=1 to stage without submitting)" >&2
    exit 1
    ;;
esac
