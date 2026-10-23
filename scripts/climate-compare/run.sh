#!/usr/bin/env bash
#
# Climate-only comparison of two par files (A and B): both run through runme with
# the ice sheet off, the forcing time evolving (tstep_method = "rel") and the
# climate updated and written every DT years. The climate 2D file of each run
# (snap.nc, esm.nc or rembo.nc on the climate grid) is then compared with
# compare.jl. See README.md.
#
# Usage (from anywhere; runs from the repo root):
#   scripts/climate-compare/run.sh [options] PAR_A PAR_B OUTDIR
#
# Options:
#   --t0 T        first forcing time [yr]          (default -130e3)
#   --t1 T        last forcing time [yr]           (default 0)
#   --pd          present day: forcing time fixed at 0 (tstep_method = "const"),
#                 two records (0 and DT), instead of --t0/--t1
#   --dt DT       climate update and output interval [yr] (default 1000)
#   -p  "k=v ..." extra runme overrides for A and B (repeatable)
#   --pa "k=v ..." extra overrides for A only (repeatable)
#   --pb "k=v ..." extra overrides for B only (repeatable)
#   --stage       stage the runs only (submit scripts written, not submitted)
#   --local       run on this machine (runme -r), no slurm
#   -e EXE        runme executable alias (default yelmox)
#
# Environment:
#   RUNOPTS       slurm options for runme (default "-q compute -w 00:30:00 --omp 16")
#
# Output: OUTDIR/A and OUTDIR/B (runme run folders), relative to the repo root.
#
# Example (GRL, paleo, ocean fraction 0.25 vs 0.5):
#   scripts/climate-compare/run.sh --pb "coef_dTo.factor=0.5" \
#       yelmox/yelmox_Greenland.nml yelmox/yelmox_Greenland.nml output/climcmp/grl
#   julia --project=scripts/climate-compare scripts/climate-compare/compare.jl output/climcmp/grl
#
set -euo pipefail

usage () { sed -n '3,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

T0=-130e3
T1=0
DT=1000
PD=0
EXE=yelmox
MODE="-rs"
P_AB=()
P_A=()
P_B=()
POS=()

while [ $# -gt 0 ]; do
    case "$1" in
        --t0)    T0=$2; shift 2 ;;
        --t1)    T1=$2; shift 2 ;;
        --dt)    DT=$2; shift 2 ;;
        --pd)    PD=1; shift ;;
        -p)      read -r -a kv <<< "$2"; P_AB+=("${kv[@]}"); shift 2 ;;
        --pa)    read -r -a kv <<< "$2"; P_A+=("${kv[@]}");  shift 2 ;;
        --pb)    read -r -a kv <<< "$2"; P_B+=("${kv[@]}");  shift 2 ;;
        --stage) MODE="-s"; shift ;;
        --local) MODE="-r"; shift ;;
        -e)      EXE=$2; shift 2 ;;
        -h|--help) usage 0 ;;
        -*)      echo "run.sh: unknown option $1" >&2; usage 1 ;;
        *)       POS+=("$1"); shift ;;
    esac
done
[ ${#POS[@]} -eq 3 ] || usage 1
PAR_A=${POS[0]}
PAR_B=${POS[1]}
OUT=${POS[2]}

cd "$(dirname "$0")/../.."     # repo root: runme needs .runme/ and the data links

RUNOPTS="${RUNOPTS:--q compute -w 00:30:00 --omp 16}"
read -r -a RUNOPTS_ARR <<< "$RUNOPTS"
[ "$MODE" = "-r" ] && RUNOPTS_ARR=()

# The group holding the timeline: [ctrl] run_step ("ctrl" = [ctrl] itself).
timeline_group () {
    local g
    g=$(awk '/^&ctrl/{f=1; next} f && /^\//{exit}
             f && $1 == "run_step" {gsub(/["'\''=]/, " "); print $2; exit}' "$1")
    echo "${g:-ctrl}"
}

# Climate-only overrides: ice sheet off, climate updated and written every DT
# (2D, 2D small and 1D output on the same cadence), no restart bundles but the
# final one.
overrides () {
    local g=$1
    local o=("comps.with_ice_sheet=False" "comps.dt_clim=$DT" "$g.dtt=$DT"
             "tm_2D.method=const" "tm_2D.dt=$DT" "tm_2Dsm.method=const" "tm_2Dsm.dt=$DT"
             "tm_1D.method=const" "tm_1D.dt=$DT" "tm_rst.method=none")
    if [ "$PD" = 1 ]; then
        o+=("$g.tstep_method=const" "$g.tstep_const=0" "$g.time_init=0" "$g.time_end=$DT")
    else
        o+=("$g.tstep_method=rel" "$g.time_init=$T0" "$g.time_end=$T1")
    fi
    echo "${o[@]}"
}

submit () {   # label par_file extra-overrides...
    local lab=$1 par=$2
    shift 2
    local g
    g=$(timeline_group "$par")
    read -r -a ov <<< "$(overrides "$g")"
    echo ">> $lab: $par -> $OUT/$lab"
    runme "$MODE" ${RUNOPTS_ARR[@]+"${RUNOPTS_ARR[@]}"} --jobname "climcmp-$lab" -e "$EXE" -n "$par" \
        -o "$OUT/$lab" -p "${ov[@]}" "$@"
}

submit A "$PAR_A" ${P_AB[@]+"${P_AB[@]}"} ${P_A[@]+"${P_A[@]}"}
submit B "$PAR_B" ${P_AB[@]+"${P_AB[@]}"} ${P_B[@]+"${P_B[@]}"}
