#!/usr/bin/env bash
#
# Nudged AWI-ESM workflow for Antarctica (yelmox, climate = esm).
#
#   Step 1  spinup     20-kyr present-day OPTIMIZED ice-sheet spin-up
#                      (sim.opt=True) -> writes a restart bundle
#   Step 2  scenarios  years 4501-5000, each branched off that bundle:
#                      ctrl          the reference climate only
#                      nudgedhigh03  plus the nudged AWI-ESM anomaly
#
# Base yelmox/yelmox_esm_Antarctica.nml with the setup nudge_ant below: ANT-16KM,
# chion PDD, and the forcing input/esm/esm_nudge.nml (Nudge_Skiba26). The
# reference climate is the mean over time_ref (1990-2019) of the RACMO2.3p2
# 1979-2022 and ORAS5 1979-2024 monthly series; experiment nudgedhigh03 adds the
# anomaly of the nudged AWI-ESM run (nudgedhigh03, the projection) minus its PI
# run (orig_PI, the ESM reference), years 3501-5000. Atmosphere only: the ocean
# has no nudged files on ANT-16KM yet.
#
# Run the steps in order on the cluster; let the spin-up finish before launching
# the scenarios (they read its restart bundle):
#
#   scripts/nudge/run_nudge_antarctica.sh spinup
#   scripts/nudge/run_nudge_antarctica.sh scenarios
#
# Stage only -- create the run dirs + SLURM submit script but do NOT submit
# (a dry run to inspect everything first): set STAGE=1
#
#   STAGE=1 scripts/nudge/run_nudge_antarctica.sh spinup
#   STAGE=1 scripts/nudge/run_nudge_antarctica.sh scenarios
#
set -euo pipefail
cd "$(dirname "$0")/../.." || exit 1               # repo root

# ---- configuration ---------------------------------------------------------
EXE="yelmox"                                       # -> libyelmox/bin/yelmox.x (climate = esm)
NML="yelmox/yelmox_esm_Antarctica.nml"
OUTROOT="output/nudge_ant"
# Scenarios: experiment | use_proj (ctrl: the reference only)
SCENARIOS=(
  "ctrl          False"
  "nudgedhigh03  True"
)

SPINUP_YEARS=20000                                 # ice-sheet opt spin-up (matches &opt cf/tf_time_end=20e3)

# The nudge setup on the base.
nudge_ant=(
    # time: the nudged years
    "transient.tstep_const=4500.0"          # [yr] time of the "const" method
    "transient.time_init=4501.0"            # [yr]
    "transient.time_end=5000.0"             # [yr]
    "transient.time_esm_ref==4501,5000"     # [yr] reference of the ESM
    "transient.time_hist==1950,1990"        # [yr] historical period
    "transient.time_proj==4501,5000"        # [yr] projection period
    # output
    "tm_2D.method=times"                    # 2D output at the start and the end
    "tm_2D.dt=5e3"
    "tm_2D.times==4501,5000"
    "tm_2Dsm.method=times"                  # small 2D output every 50 yr
    "tm_2Dsm.dt=1.0"
    "tm_2Dsm.times==4501,4550,4600,4650,4700,4750,4800,4850,4900,4950,5000"
    # grid
    "domain.grid_hub=ANT-16KM"              # hub, and through it every component grid
    # relaxation and optimization
    "relax.tau2=6000.0"                     # [yr] relaxation timescale at time2
    "relax.time1=1e3"                       # [yr] end of the tau1 period
    "relax.time2=1e3"                       # [yr] end of the relaxation
    "opt.tf_min=-3.0"                       # [degC] minimum tf_corr
    # surface: chion PDD
    "chion.model=pdd"                       # snowpack model
    "chion.nml_pdd=pdd_chion"               # its parameters: the base's [pdd_chion] (cg05)
    "smbpal.abl_method=pdd"                 # smbpal ablation
    # ocean: PICO
    "pico.gamma_tstar=630.0"                # [m/yr] turbulent temperature exchange velocity (Reese et al., 2018)
    "pico.C_over=1.0"                       # [Sv m3/kg] overturning strength (Reese et al., 2018)
    # ESM forcing
    "esm.par_file=input/esm/esm_nudge.nml"
    "esm.esm_name=AWI-ESM"                  # {gcm}; not in the nudge file names
    "esm.use_esm=True"                      # ESM snapshots
)

# runme submit options. STAGE=1 writes the submit script without submitting.
if [ "${STAGE:-0}" = 1 ]; then SUBMIT="-s"; else SUBMIT="-rs"; fi
# HPCOPT_SPINUP / HPCOPT_SCEN can be set in the environment (other queues, e.g. on
# albedo: HPCOPT_SPINUP="-q 48h -w 48:00:00 --omp 8").
HPCOPT_SPINUP="${HPCOPT_SPINUP:--q shared -w 7-00:00:00 -m 20G --omp 16}"
HPCOPT_SCEN="${HPCOPT_SCEN:--q compute -w 08:00:00 --omp 16}"

SPINUP_OUT="$OUTROOT/spinup"
# The spin-up's final restart bundle. yelmox names it restart-<time/1e3 %.3f>-kyr.
# Absolute path: the executable runs from inside the scenario's run dir, so a
# repo-root-relative path would not resolve.
BUNDLE="$(pwd)/$SPINUP_OUT/restart-$(awk "BEGIN{printf \"%.3f\", $SPINUP_YEARS/1000}")-kyr"

# ---- steps -----------------------------------------------------------------
case "${1:-}" in
  spinup)
    runme $SUBMIT $HPCOPT_SPINUP -e "$EXE" -n "$NML" -o "$SPINUP_OUT" \
      -p "${nudge_ant[@]}" \
         ctrl.run_step=spinup \
         spinup.time_init=0 spinup.time_end="$SPINUP_YEARS"
    ;;
  scenarios)
    for row in "${SCENARIOS[@]}"; do
      read -r exp proj <<<"$row"
      runme $SUBMIT $HPCOPT_SCEN -e "$EXE" -n "$NML" -o "$OUTROOT/$exp" \
        -p "${nudge_ant[@]}" \
           ctrl.run_step=transient esm.experiment="$exp" esm.use_proj="$proj" \
           ctrl.restart="$BUNDLE"
    done
    ;;
  *)
    echo "usage: $0 {spinup|scenarios}   (prefix STAGE=1 to stage without submitting)" >&2
    exit 1
    ;;
esac
