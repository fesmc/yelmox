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
# SCAFFOLD STATUS: the RAMP (esm-up2p0, 232 yr) and BOTH STABILISATIONS
# (esm-up2p0-gwl2p0 / -gwl4p0, 50 yr) are wired end-to-end for IPSL-CM6-ESMCO2 @
# GRL-8KM. The stabilisations start already warm, so they use a separate par_file
# (esm_grl_tipmip_gwl.nml) that pins esm_ref to the ramp's piControl branch point
# (zero reference) and uses a 50-yr projection clock; the scenarios step below
# overrides esm.par_file + transient.time_end per experiment (see the RUNS table).
# NOTE: the second forcing model EC-Earth3-ESM-1 (GRL-4KM only) is NOT wired here --
# GRL-4KM currently lacks the ORAS4 ocean reference (empty dir) and a combined-pr
# MAR climatology on this machine, so it cannot run the forcing setup as-is.
#
set -euo pipefail
cd "$(dirname "$0")/.." || exit 1                  # repo root

# ---- configuration ---------------------------------------------------------
EXE="esm"                                          # -> libyelmox/bin/yelmox_esm.x
NML="yelmox_esm/yelmox_esm_Greenland_tipmip.nml"
OUTROOT="output/tipmip_grl"
GCM="IPSL-CM6-ESMCO2"
GRID="GRL-8KM"
SPINUP_YEARS=10                                    # forcing-only: short suffices

# Scenario table:  experiment | par_file (-> esm.par_file) | transient.time_end
# Ramp uses the base par nml (232 yr); stabilisations use the gwl par nml (50 yr,
# ramp-referenced zero baseline). See the esm_grl_tipmip*.nml headers.
RUNS=(
  "esm-up2p0         input/esm/esm_grl_tipmip.nml      232"
  "esm-up2p0-gwl2p0  input/esm/esm_grl_tipmip_gwl.nml   50"
  "esm-up2p0-gwl4p0  input/esm/esm_grl_tipmip_gwl.nml   50"
)

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
    # Ramp then both stabilisations, each restart-branched off the spin-up bundle.
    for row in "${RUNS[@]}"; do
      read -r exp par tend <<<"$row"
      runme $SUBMIT $HPCOPT -e "$EXE" -n "$NML" -o "$OUTROOT/$exp" \
        -p ctrl.run_step=transient esm.experiment="$exp" esm.esm_name="$GCM" \
           esm.par_file="$par" esm.use_esm=True esm.use_hist=False esm.use_proj=True \
           yelmo.grid_name="$GRID" htopo.grid_name="$GRID" \
           transient.time_end="$tend" \
           coupling.restart="$BUNDLE"
    done
    ;;
  *)
    echo "usage: $0 {spinup|scenarios}   (prefix STAGE=1 to stage without submitting)" >&2
    exit 1
    ;;
esac
