#!/usr/bin/env bash
#
# Smoke tests for the yelmox_bipolar driver: the coupling setups of
# yelmox_bipolar/params.sh, each a short run (the par file's time_init to
# time_end) that must reach "yelmox_bipolar: run complete".
#
# Runs LOCALLY through the usual runme workflow: `-r` runs the executable
# directly, one setup after the other. Each setup writes to its own directory
# under $OUT. Run by hand with:
#
#     bash yelmox_bipolar/run_tests.sh                     # every setup (bip_setups)
#     bash yelmox_bipolar/run_tests.sh bip_full bip_gis_obm
#     bash yelmox_bipolar/run_tests.sh bip_full+bip_ocean_anomaly   # a setup + a modifier
#
# Prerequisites:
#   - ice-sheet input data for the north (Greenland, GRL-16KM) and the south
#     (Antarctica, ANT-32KM) under ice_data/ and fesmdata/.
#   - for the setups with couple_fwf_north/south: the hydrographic masks of
#     [ctrl] hydro_mask_north/south (…_HYDRO-BASINS_zickfeld2004.nc).

cd "$(dirname "$0")/.." || exit 1        # repo root
source yelmox_bipolar/params.sh

NML="yelmox_bipolar/yelmox_bipolar_Bipolar.nml"
EXE="bipolar"
OUT="${OUT:-tmp/bipolar-tests}"
TIMEOUT="${TIMEOUT:-3600}"               # per-setup wall-clock cap [s]

setups=("$@")
[ ${#setups[@]} -eq 0 ] && setups=("${bip_setups[@]}" "bip_full+bip_ocean_anomaly")

echo "=== building yelmox_bipolar ==="
make yelmox_bipolar || { echo "build failed"; exit 1; }

status=0
summary=()
run_test() {
    # run_test <setup>: the -p arrays named in <setup>, joined by "+"
    local setup="$1" outdir="$OUT/${1//+/_}" name
    local -a params=()
    IFS="+" read -r -a names <<< "$setup"
    for name in "${names[@]}"; do
        if ! declare -p "$name" >/dev/null 2>&1; then
            echo "  [FAIL] $setup (no array $name in yelmox_bipolar/params.sh)"
            summary+=("FAIL $setup"); status=1; return
        fi
        eval 'params+=("${'"$name"'[@]}")'
    done

    echo
    echo "=================================================================="
    echo "  $setup"
    echo "  -> $outdir"
    echo "  -p ${params[*]}"
    echo "=================================================================="
    rm -rf "$outdir"                     # fresh dir so the test is re-runnable

    # runme -r launches the executable in the BACKGROUND (…/out.out &) and returns
    # immediately, so poll the run log for a clean finish vs. an abort. Keep this
    # shell alive meanwhile, otherwise the backgrounded job is killed on exit.
    local -a pargs=()
    [ ${#params[@]} -gt 0 ] && pargs=(-p "${params[@]}")
    if ! runme -r -e "$EXE" -n "$NML" -o "$outdir" "${pargs[@]}"; then
        echo "  [FAIL] $setup (runme could not stage/launch the run)"
        summary+=("FAIL $setup"); status=1; return
    fi

    local log="$outdir/out.out" waited=0
    while true; do
        if grep -q "yelmox_bipolar: run complete" "$log" 2>/dev/null; then
            echo "  [PASS] $setup"; summary+=("PASS $setup"); return
        fi
        if grep -qE "STOP [0-9]+|Fortran runtime error|Backtrace for|forrtl: severe" "$log" 2>/dev/null; then
            echo "  [FAIL] $setup (aborted, see $log)"
            grep -m3 -E ":: error|STOP [0-9]+|forrtl" "$log" | sed 's/^/         /'
            summary+=("FAIL $setup"); status=1; return
        fi
        if ! pgrep -f "yelmox_bipolar.x" >/dev/null 2>&1; then
            echo "  [FAIL] $setup (process exited without a completion marker; see $log)"
            summary+=("FAIL $setup"); status=1; return
        fi
        sleep 10; waited=$((waited + 10))
        if [ "$waited" -ge "$TIMEOUT" ]; then
            echo "  [FAIL] $setup (timeout after ${TIMEOUT}s; still running, see $log)"
            summary+=("FAIL $setup"); status=1; return
        fi
    done
}

for setup in "${setups[@]}"; do
    run_test "$setup"
done

echo
echo "=== summary ==="
printf '  %s\n' "${summary[@]}"
exit "$status"
