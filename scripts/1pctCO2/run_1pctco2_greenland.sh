#!/usr/bin/env bash
#
# 1pctCO2 workflow for Greenland (yelmox, climate = esm).
#
#   Step 1  spinup     15-kyr present-day OPTIMIZED ice-sheet spin-up
#                      (sim.opt=True, as ISMIP7) -> writes a restart bundle
#   Step 2  scenarios  the 1pctCO2 run, branched off that bundle
#
# The ice sheet + isostasy are ACTIVE. 1pctCO2 = idealized CMIP experiment,
# atmospheric CO2 +1%/yr to 4xCO2 at ~yr 140. The forcing is the GCM's ABSOLUTE
# tas/pr/thetao/so, self-referenced to the run start (~1xCO2), read from the global
# CMIP files and remapped online onto the run's grid; see input/esm/esm_grl_1pctCO2.nml.
# The reference climatology is MAR (atmosphere) + ERA-INT-ORAS4 (ocean).
#
# Run the steps in order on the cluster; let the spin-up finish first:
#
#   scripts/1pctCO2/run_1pctco2_greenland.sh spinup
#   scripts/1pctCO2/run_1pctco2_greenland.sh scenarios
#
# The GCM is set with the GCM env var (default MPI-ESM1-2-LR; the spin-up does not
# depend on it): MPI-ESM1-2-LR | NorESM2-MM | UKESM1-0-LL | IPSL-CM6A-LR | CESM2
#
# Stage only (create dirs + SLURM submit script, do NOT submit): STAGE=1
#
#   STAGE=1 scripts/1pctCO2/run_1pctco2_greenland.sh spinup
#
set -euo pipefail
cd "$(dirname "$0")/../.." || exit 1               # repo root

# ---- configuration ---------------------------------------------------------
EXE="yelmox"                                       # -> libyelmox/bin/yelmox.x (climate = esm)
NML="yelmox/yelmox_esm_Greenland_1pctCO2.nml"
OUTROOT="output/1pctco2_grl"
GCM="${GCM:-MPI-ESM1-2-LR}"
GRID="GRL-8KM"

SPINUP_YEARS=15000                                 # ice-sheet opt spin-up (matches &opt cf/tf_time_end=15e3)
PROJ_INIT=2020                                     # 1pctCO2 start year (year 0 = 1xCO2)
PROJ_END=2160                                      # 1pctCO2 end year (140 yr)

# runme submit options. STAGE=1 writes the submit script without submitting.
if [ "${STAGE:-0}" = 1 ]; then SUBMIT="-s"; else SUBMIT="-rs"; fi
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
      -p ctrl.run_step=spinup sim.relax=True sim.opt=True esm.experiment=ctrl \
         domain.grid_hub="$GRID" \
         spinup.time_init=0 spinup.time_end="$SPINUP_YEARS"
    ;;
  scenarios)
    runme $SUBMIT $HPCOPT_SCEN -e "$EXE" -n "$NML" -o "$OUTROOT/1pctCO2_$GCM" \
      -p ctrl.run_step=transient esm.experiment=1pctCO2 esm.esm_name="$GCM" \
         esm.use_esm=True esm.use_hist=False esm.use_proj=True \
         domain.grid_hub="$GRID" \
         ctrl.restart="$BUNDLE" \
         transient.time_init="$PROJ_INIT" transient.time_end="$PROJ_END"
    ;;
  *)
    echo "usage: [GCM=...] $0 {spinup|scenarios}   (prefix STAGE=1 to stage without submitting)" >&2
    exit 1
    ;;
esac
