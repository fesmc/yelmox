#!/usr/bin/env bash
#
# TIPMIP forcing-only workflow for Greenland (refactored yelmox_esm) -- SCAFFOLD.
#
#   Step 1  spinup     reference-climate spin-up  -> writes a restart bundle
#   Step 2  scenarios  TIPMIP experiment(s), each branched off that bundle
#
# All runs are forcing-only (coupling.with_ice_sheet=False, with_isostasy=False):
# they produce climate + ocean forcing from the TIPMIP piControl anomalies
# (tas_anomaly / pr_ratio / TF_anomaly). See input/esm/esm_grl_tipmip.nml and its
# header for the anomaly-referencing assumptions and scaffold caveats.
#
# Run the steps in order on the cluster; let the spin-up finish first:
#
#   yelmox_esm/run_tipmip_greenland.sh spinup
#   yelmox_esm/run_tipmip_greenland.sh scenarios
#
# Stage only (create dirs + SLURM submit script, do NOT submit): STAGE=1
#
#   STAGE=1 yelmox_esm/run_tipmip_greenland.sh spinup
#
# SCAFFOLD STATUS: only the RAMP (esm-up2p0) is wired end-to-end -- its 0..232 yr
# clock and start-referenced anomaly match the run config. The stabilisations
# (esm-up2p0-gwl2p0 / -gwl4p0) are listed but commented: they span different
# piControl years (109..159 / 232..282) and start already warm, so they need the
# &transient time_* and a proper zero reference set before enabling (see the
# esm_grl_tipmip.nml header). Second forcing model EC-Earth3-ESM-1 exists only at
# GRL-4KM (set GCM + GRID accordingly).
#
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1                  # repo root

# ---- configuration ---------------------------------------------------------
EXE="esm"                                          # -> libyelmox/bin/yelmox_esm.x
NML="yelmox_esm/yelmox_esm_Greenland_tipmip.nml"
OUTROOT="output/tipmip_grl"
GCM="IPSL-CM6-ESMCO2"
GRID="GRL-8KM"
SCENARIOS=(esm-up2p0)                              # ramp; add stabilisations once tuned:
                                                   #   esm-up2p0-gwl2p0 esm-up2p0-gwl4p0
SPINUP_YEARS=10                                    # forcing-only: short suffices

# runme submit options. STAGE=1 writes the submit script without submitting.
if [ "${STAGE:-0}" = 1 ]; then SUBMIT="-s"; else SUBMIT="-rs"; fi
HPCOPT="-q compute -w 02:00:00 --omp 8"

SPINUP_OUT="$OUTROOT/spinup"
# Absolute path: the executable runs from inside the scenario's run dir.
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
           coupling.restart="$BUNDLE"
    done
    ;;
  *)
    echo "usage: $0 {spinup|scenarios}   (prefix STAGE=1 to stage without submitting)" >&2
    exit 1
    ;;
esac
