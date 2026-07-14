#!/usr/bin/env bash
#
# ISMIP7 forcing-only workflow for Antarctica (refactored yelmox_esm).
#
#   Step 1  spinup     reference-climate spin-up  -> writes a restart bundle
#   Step 2  scenarios  ssp585, branched off that bundle
#
# All runs are forcing-only (coupling.with_ice_sheet=False, with_isostasy=False):
# they produce ISMIP7 climate + ocean forcing, no ice dynamics.
#
# Run the steps in order on the cluster; let the spin-up finish before launching
# the scenarios (they read its restart bundle):
#
#   yelmox_esm/run_ismip7_antarctica.sh spinup
#   yelmox_esm/run_ismip7_antarctica.sh scenarios
#
# Stage only -- create the run dirs + SLURM submit script but do NOT submit
# (a dry run to inspect everything first): set STAGE=1
#
#   STAGE=1 yelmox_esm/run_ismip7_antarctica.sh spinup
#   STAGE=1 yelmox_esm/run_ismip7_antarctica.sh scenarios
#
# NOTE: the ANT-8KM grid needs a large stack. The SLURM submit script sets it,
# but for a local run first do:  ulimit -s unlimited
#
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1                  # repo root

# ---- configuration ---------------------------------------------------------
EXE="esm"                                          # -> libyelmox/bin/yelmox_esm.x
NML="yelmox_esm/yelmox_esm_Antarctica_ismip7.nml"
OUTROOT="output/ismip7_ant"
GCM="CESM2-WACCM"
GRID="ANT-8KM"
SCENARIOS=(ssp585)                                 # only ssp585 present on Levante

SPINUP_YEARS=10                                    # forcing-only: short suffices
PROJ_END=2300                                      # scenario end year (CE)

# runme submit options. STAGE=1 writes the submit script without submitting.
if [ "${STAGE:-0}" = 1 ]; then SUBMIT="-s"; else SUBMIT="-rs"; fi
HPCOPT="-q compute -w 08:00:00 --omp 16"           # ANT-8KM is heavy; give it threads

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
