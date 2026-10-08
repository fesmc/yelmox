#!/usr/bin/env bash
#
# TIPMIP workflow for Greenland (yelmox, climate = esm).
#
#   Step 1  spinup          15-kyr present-day OPTIMIZED ice-sheet spin-up
#                           (sim.opt=True, as ISMIP7) -> writes a restart bundle
#   Step 2  ramp            esm-up2p0 (232 yr), branched off the spin-up bundle;
#                           also writes restart bundles at the stabilisation
#                           branch years
#   Step 3  stabilisations  esm-up2p0-gwl2p0 / -gwl4p0 (50 yr), each branched off
#                           the ramp's bundle at its branch year
#
# The ice sheet + isostasy are ACTIVE. The forcing is the TIPMIP piControl anomalies
# (tas_anomaly / pr_ratio / TF_anomaly) on top of the MAR + ERA-INT-ORAS4 reference
# climatology; TIPMIP has no salinity (dso = 0) or subglacial discharge (Qd = 0).
# See input/esm/esm_grl_tipmip.nml (ramp) and esm_grl_tipmip_gwl.nml (stabilisations).
#
# Run the steps in order on the cluster; let each finish before the next (they read
# the previous step's restart bundles):
#
#   scripts/tipmip/run_tipmip_greenland.sh spinup
#   scripts/tipmip/run_tipmip_greenland.sh ramp
#   scripts/tipmip/run_tipmip_greenland.sh stabilisations
#
# Stage only (create dirs + SLURM submit script, do NOT submit): STAGE=1
#
#   STAGE=1 scripts/tipmip/run_tipmip_greenland.sh spinup
#
# MODELS (select with the MODEL env var, default ipsl):
#   MODEL=ipsl     IPSL-CM6-ESMCO2 @ GRL-8KM -- ramp + both stabilisations
#   MODEL=ecearth  EC-Earth3-ESM-1 @ GRL-4KM -- ramp only (no stabilisations)
#
# The stabilisations start already warm, so they use a separate par_file
# (esm_grl_tipmip_gwl.nml) that pins esm_ref to the ramp's piControl branch point
# (zero reference) with a 50-yr projection clock (years 0..49). They branch off the
# ramp at the years their forcing continues from: gwl2p0 at ramp year 109,
# gwl4p0 at ramp year 232 (the ramp's end).
#
set -euo pipefail
cd "$(dirname "$0")/../.." || exit 1               # repo root

# ---- configuration ---------------------------------------------------------
EXE="yelmox"                                       # -> libyelmox/bin/yelmox.x (climate = esm)
NML="yelmox/yelmox_esm_Greenland_tipmip.nml"
SPINUP_YEARS=15000                                 # ice-sheet opt spin-up (matches &opt cf/tf_time_end=15e3)
RAMP_YEARS=232                                     # esm-up2p0
STAB_YEARS=50                                      # esm-up2p0-gwl*
MODEL="${MODEL:-ipsl}"                             # ipsl | ecearth

# Per-model config + stabilisation table (experiment | ramp branch year).
case "$MODEL" in
  ipsl)
    GCM="IPSL-CM6-ESMCO2"; GRID="GRL-8KM"; OUTROOT="output/tipmip_grl_ipsl"
    HPCOPT_SPINUP="-q compute -w 08:00:00 --omp 8"
    HPCOPT_SCEN="-q compute -w 02:00:00 --omp 8"
    STABS=(
      "esm-up2p0-gwl2p0  109"
      "esm-up2p0-gwl4p0  232"
    )
    ;;
  ecearth)
    GCM="EC-Earth3-ESM-1"; GRID="GRL-4KM"; OUTROOT="output/tipmip_grl_ecearth"
    HPCOPT_SPINUP="-q compute -w 08:00:00 --omp 32"   # GRL-4KM: verify the spin-up fits the queue limit
    HPCOPT_SCEN="-q compute -w 08:00:00 --omp 16"
    STABS=()
    ;;
  *) echo "unknown MODEL='$MODEL' (use: ipsl | ecearth)" >&2; exit 1 ;;
esac

# runme submit options. STAGE=1 writes the submit script without submitting.
if [ "${STAGE:-0}" = 1 ]; then SUBMIT="-s"; else SUBMIT="-rs"; fi

# yelmox names restart bundles restart-<time/1e3 %.3f>-kyr. Absolute paths: the
# executable runs from inside the run dir, so a repo-root-relative path would not resolve.
bundle() { echo "$(pwd)/$1/restart-$(awk "BEGIN{printf \"%.3f\", $2/1000}")-kyr"; }
SPINUP_OUT="$OUTROOT/spinup"
RAMP_OUT="$OUTROOT/esm-up2p0"

# Ramp restart times: every stabilisation branch year (the end is always written).
RST_TIMES=""
for row in ${STABS[@]+"${STABS[@]}"}; do
  read -r _ yr <<<"$row"
  RST_TIMES="${RST_TIMES:+$RST_TIMES,}$yr"
done

# ---- steps -----------------------------------------------------------------
case "${1:-}" in
  spinup)
    runme $SUBMIT $HPCOPT_SPINUP -e "$EXE" -n "$NML" -o "$SPINUP_OUT" \
      -p ctrl.run_step=spinup sim.relax=True sim.opt=True \
         esm.experiment=ctrl esm.esm_name="$GCM" \
         domain.grid_hub="$GRID" \
         spinup.time_init=0 spinup.time_end="$SPINUP_YEARS"
    ;;
  ramp)
    rst=()
    # "==" passes the list as one vector value (a plain comma list is a runme ensemble)
    [ -n "$RST_TIMES" ] && rst=(tm_rst.method=times "tm_rst.times==$RST_TIMES")
    runme $SUBMIT $HPCOPT_SCEN -e "$EXE" -n "$NML" -o "$RAMP_OUT" \
      -p ctrl.run_step=transient esm.experiment=esm-up2p0 esm.esm_name="$GCM" \
         esm.par_file=input/esm/esm_grl_tipmip.nml \
         esm.use_esm=True esm.use_hist=False esm.use_proj=True \
         domain.grid_hub="$GRID" \
         transient.time_init=0 transient.time_end="$RAMP_YEARS" \
         ctrl.restart="$(bundle "$SPINUP_OUT" "$SPINUP_YEARS")" ${rst[@]+"${rst[@]}"}
    ;;
  stabilisations)
    if [ ${#STABS[@]} -eq 0 ]; then echo "MODEL=$MODEL has no stabilisations." >&2; exit 1; fi
    for row in "${STABS[@]}"; do
      read -r exp yr <<<"$row"
      runme $SUBMIT $HPCOPT_SCEN -e "$EXE" -n "$NML" -o "$OUTROOT/$exp" \
        -p ctrl.run_step=transient esm.experiment="$exp" esm.esm_name="$GCM" \
           esm.par_file=input/esm/esm_grl_tipmip_gwl.nml \
           esm.use_esm=True esm.use_hist=False esm.use_proj=True \
           domain.grid_hub="$GRID" \
           transient.time_init=0 transient.time_end="$STAB_YEARS" \
           ctrl.restart="$(bundle "$RAMP_OUT" "$yr")"
    done
    ;;
  *)
    echo "usage: [MODEL=ipsl|ecearth] $0 {spinup|ramp|stabilisations}   (prefix STAGE=1 to stage without submitting)" >&2
    exit 1
    ;;
esac
