#!/usr/bin/env bash
#
# 1pctCO2 forcing-only workflow for Antarctica (refactored yelmox_esm) -- SCAFFOLD.
#
#   Step 1  spinup     reference-climate spin-up  -> writes a restart bundle
#   Step 2  scenarios  the 1pctCO2 run, branched off that bundle
#
# All runs are forcing-only (coupling.with_ice_sheet=False, with_isostasy=False):
# they produce climate + ocean forcing from the 1pctCO2 CMIP fields, no ice dynamics.
# 1pctCO2 = idealized CMIP experiment, atmospheric CO2 +1%/yr to 4xCO2 at ~yr 140.
# Forcing is ABSOLUTE tas/pr/thetao/so, self-referenced to the run start (~1xCO2);
# see input/esm/esm_ant_1pctCO2.nml and its header.
#
# *** UNTESTED / DATA-PENDING ***  The 1pctCO2 forcing is NOT on this machine yet
# (the whole ice_data/1pctCO2/ tree is absent). This scaffold is wired against the
# expected layout (realization r1i1p4f1, window 2020-2160) so it should run once the
# data is staged. Until then only STAGE=1 (dir/submit-script staging) will succeed;
# a real run will stop in varslice when it cannot find the 1pctCO2 files.
#
# Run the steps in order on the cluster; let the spin-up finish first:
#
#   yelmox_esm/run_1pctco2_antarctica.sh spinup
#   yelmox_esm/run_1pctco2_antarctica.sh scenarios
#
# Stage only (create dirs + SLURM submit script, do NOT submit): STAGE=1
#
#   STAGE=1 yelmox_esm/run_1pctco2_antarctica.sh spinup
#
# NOTE: larger ANT grids need a big stack. The SLURM submit script sets it; for a
# local run first do:  ulimit -s unlimited
#
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1                  # repo root

# ---- configuration ---------------------------------------------------------
EXE="esm"                                          # -> libyelmox/bin/yelmox_esm.x
NML="yelmox_esm/yelmox_esm_Antarctica_1pctCO2.nml"
OUTROOT="output/1pctco2_ant"
GCM="1pctCO2-r1i1p4f1"                             # label only (no {gcm} templating in the par nml)
GRID="ANT-32KM"                                    # lightest grid; bump once data + grid confirmed
SCENARIOS=(1pctCO2)

SPINUP_YEARS=10                                    # forcing-only: short suffices
PROJ_INIT=2020                                     # 1pctCO2 start year (run year 0 = 1xCO2)
PROJ_END=2160                                      # 1pctCO2 end year (140 yr)

# runme submit options. STAGE=1 writes the submit script without submitting.
if [ "${STAGE:-0}" = 1 ]; then SUBMIT="-s"; else SUBMIT="-rs"; fi
HPCOPT="-q compute -w 04:00:00 --omp 16"

SPINUP_OUT="$OUTROOT/spinup"
# Absolute path: the executable runs from inside the scenario's run dir, so a
# repo-root-relative path would not resolve.
BUNDLE="$(pwd)/$SPINUP_OUT/restart-$(awk "BEGIN{printf \"%.3f\", $SPINUP_YEARS/1000}")-kyr"

# ---- steps -----------------------------------------------------------------
case "${1:-}" in
  spinup)
    runme $SUBMIT $HPCOPT -e "$EXE" -n "$NML" -o "$SPINUP_OUT" \
      -p ctrl.run_step=spinup esm.experiment=ctrl esm.esm_name="$GCM" \
         yelmo.grid_name="$GRID" htopo.grid_name="$GRID" \
         spinup.time_init=0 spinup.time_end="$SPINUP_YEARS"
    ;;
  scenarios)
    for exp in "${SCENARIOS[@]}"; do
      runme $SUBMIT $HPCOPT -e "$EXE" -n "$NML" -o "$OUTROOT/$exp" \
        -p ctrl.run_step=transient esm.experiment="$exp" esm.esm_name="$GCM" \
           esm.use_esm=True esm.use_hist=False esm.use_proj=True \
           yelmo.grid_name="$GRID" htopo.grid_name="$GRID" \
           coupling.restart="$BUNDLE" \
           transient.time_init="$PROJ_INIT" transient.time_end="$PROJ_END"
    done
    ;;
  *)
    echo "usage: $0 {spinup|scenarios}   (prefix STAGE=1 to stage without submitting)" >&2
    exit 1
    ;;
esac
