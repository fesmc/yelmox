#!/usr/bin/env bash
#
# A long run as a chain of short jobs ("chunks"), for queues where short jobs
# start sooner: chunk k runs from T_k to T_k+1 in its own run folder
# OUTDIR/cNN and restarts from the bundle chunk k-1 wrote at its time_end
# (restart-<T_k/1000>-kyr); each job waits for the previous one (Slurm
# afterok). All chunks are staged and submitted at once. See docs/run-chunks.md.
#
# Usage (from anywhere; runs from the repo root):
#   scripts/run_chunks.sh [options] PAR OUTDIR T0 T1 [T2 ...]
#
#   PAR      par file
#   OUTDIR   output folder; chunk k in OUTDIR/cNN
#   T0 T1 .. chunk boundaries [yr] in the time of the timeline group
#            (time_init of chunk k = T_k, time_end = T_k+1)
#
# Options:
#   -r BUNDLE   restart bundle of the first chunk (default: cold start or the
#               par file's ctrl.restart)
#   -d JOBID    the first chunk waits for this job (afterok)
#   -g GROUP    timeline group (default: [ctrl] run_step of PAR, else ctrl)
#   -p "k=v .." overrides for all chunks (repeatable)
#   --p-first "k=v .."  overrides for the first chunk only (repeatable)
#   --p-rest  "k=v .."  overrides for the later chunks only (repeatable), e.g.
#               sim.relax=False ytrc.elsa_restart=True
#   -e EXE      runme executable alias (default yelmox)
#   -n NAME     job name prefix (default: the name of OUTDIR)
#
# Environment:
#   RUNOPTS     slurm options for runme (default "-q compute -w 00:30:00 --omp 16")
#
# Output: the job of each chunk on stderr; on stdout the job id of the last
# chunk and its final bundle ("JOBID BUNDLE"), to chain a next stage with -d/-r.
#
set -euo pipefail

usage () { sed -n '3,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

RST=""; DEP=""; G=""; EXE=yelmox; NAME=""
P_ALL=(); P_FIRST=(); P_REST=(); POS=()
while [ $# -gt 0 ]; do
    case "$1" in
        -r) RST=$2; shift 2 ;;
        -d) DEP=$2; shift 2 ;;
        -g) G=$2; shift 2 ;;
        -e) EXE=$2; shift 2 ;;
        -n) NAME=$2; shift 2 ;;
        -p)        read -r -a kv <<< "$2"; P_ALL+=("${kv[@]}");   shift 2 ;;
        --p-first) read -r -a kv <<< "$2"; P_FIRST+=("${kv[@]}"); shift 2 ;;
        --p-rest)  read -r -a kv <<< "$2"; P_REST+=("${kv[@]}");  shift 2 ;;
        -h|--help) usage 0 ;;
        -[0-9]*) POS+=("$1"); shift ;;       # negative chunk boundary
        -*) echo "run_chunks.sh: unknown option $1" >&2; usage 1 ;;
        *)  POS+=("$1"); shift ;;
    esac
done
[ ${#POS[@]} -ge 4 ] || usage 1
PAR=${POS[0]}; OUT=${POS[1]}; TIMES=("${POS[@]:2}")

cd "$(dirname "$0")/.."        # repo root: runme needs .runme/ and the data links

RUNOPTS="${RUNOPTS:--q compute -w 00:30:00 --omp 16}"
read -r -a RUNOPTS_ARR <<< "$RUNOPTS"
[ -n "$NAME" ] || NAME=$(basename "$OUT")
if [ -z "$G" ]; then
    G=$(awk '/^&ctrl/{f=1; next} f && /^\//{exit}
             f && $1 == "run_step" {gsub(/["'\''=]/, " "); print $2; exit}' "$PAR")
    G=${G:-ctrl}
fi
OUT=$(realpath -m "$OUT")      # restart paths must be absolute (runs start in their folder)

# The bundle a run writes at time t (restart_bundle_dir of libs/kryos_startup.f90).
bundle () { printf "restart-%.3f-kyr" "$(echo "$1/1000" | bc -l)"; }

jid=""
for ((k=0; k<${#TIMES[@]}-1; k++)); do
    t0=${TIMES[k]}; t1=${TIMES[k+1]}
    dir=$(printf "%s/c%02d" "$OUT" "$k")
    p=("$G.time_init=$t0" "$G.time_end=$t1" ${P_ALL[@]+"${P_ALL[@]}"})
    [ -n "$RST" ] && p+=("ctrl.restart=$RST")
    if [ $k -eq 0 ]; then p+=(${P_FIRST[@]+"${P_FIRST[@]}"}); else p+=(${P_REST[@]+"${P_REST[@]}"}); fi
    runme -s "${RUNOPTS_ARR[@]}" --jobname "$NAME$k" -e "$EXE" -n "$PAR" -o "$dir" -p "${p[@]}" > /dev/null
    if [ -n "$DEP" ]; then
        jid=$(cd "$dir" && sbatch --parsable --dependency=afterok:"$DEP" job.submit)
    else
        jid=$(cd "$dir" && sbatch --parsable job.submit)
    fi
    echo "$NAME c$k [$t0, $t1]: job $jid" >&2
    DEP=$jid
    RST="$dir/$(bundle "$t1")"
done
echo "$jid $RST"
