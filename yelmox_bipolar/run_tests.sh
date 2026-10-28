#!/usr/bin/env bash
#
# Smoke tests for the yelmox_bipolar driver: the coupling setups of
# yelmox_bipolar/params.sh, each a short run (the par file's time_init to
# time_end) that must reach "yelmox_bipolar: run complete", and the
# consistency checks at init, each a setup that must stop with its error.
# The forcing tests (argument `forcing`): the setups of bip_forcing_setups,
# ice sheets frozen and gmt a ramp, each checked by
# yelmox_bipolar/check_forcing.jl (the coefficients in forcing_ts.nc, the
# climate of each domain, the box model), the last one also restarted from
# its bundle at 50 years and compared with its continuous run.
#
# Runs LOCALLY through the usual runme workflow: `-r` runs the executable
# directly, one test after the other. Each test writes to its own directory
# under $OUT. Run by hand with:
#
#     bash yelmox_bipolar/run_tests.sh                     # every setup and check
#     bash yelmox_bipolar/run_tests.sh bip_full bip_gis_obm
#     bash yelmox_bipolar/run_tests.sh bip_full+bip_ocean_anomaly   # a setup + a modifier
#     bash yelmox_bipolar/run_tests.sh forcing             # the forcing tests
#
# Prerequisites:
#   - ice-sheet input data for the north (Greenland, GRL-16KM) and the south
#     (Antarctica, ANT-32KM) under ice_data/ and fesmdata/.
#   - for the setups with couple_fwf_north/south: the hydrographic masks of
#     [ctrl] hydro_mask_north/south (…_HYDRO-BASINS_zickfeld2004.nc).
#   - for the forcing tests: Julia with the environment of yelmox_bipolar
#     (NCDatasets; $JULIA, default julia).

cd "$(dirname "$0")/.." || exit 1        # repo root
source yelmox_bipolar/params.sh

NML="yelmox_bipolar/yelmox_bipolar_Bipolar.nml"
EXE="bipolar"
OUT="${OUT:-tmp/bipolar-tests}"
TIMEOUT="${TIMEOUT:-3600}"               # per-test wall-clock cap [s]
JULIA="${JULIA:-julia}"

# The consistency checks of obm_coupling_init: overrides that break a rule,
# after a setup, and the error each must stop with.
chk_other_box=("forcing_north.terms==regional_north,obm_south")
chk_box_and_air=("forcing.coefs==gmt,dTa_north,dTa_tropics,dTa_south,dTo,obm_north,obm_south"
                 "forcing_north.terms==regional_north,ocean,obm_north")
chk_shared_air=("forcing_south.terms==regional_south,ocean")
chk_box_unforced=("ctrl.active_obm=True" "coef_gmt.method=ramp-slope")
chk_unused=("forcing.coefs==gmt,dTa_north,dTa_tropics,dTa_south,dTo")
checks=("bip_full+chk_other_box"         "is derived from the box obm_south"
        "bip_full+chk_box_and_air"       "(its ocean from the box) and"
        "bip_uncoupled+chk_shared_air"   "is used by the terms of both domains"
        "bip_uncoupled+chk_box_unforced" "box model is not forced by the atmosphere"
        "bip_uncoupled+chk_unused"       "'dTa_tropics' is in [forcing] coefs, but no term")

setups=("$@")
run_checks=0
run_forcing=0
if [ ${#setups[@]} -eq 0 ]; then
    setups=("${bip_setups[@]}" "bip_full+bip_ocean_anomaly")
    run_checks=1
elif [ "${setups[0]}" = "forcing" ]; then
    setups=("${bip_forcing_setups[@]}")
    run_forcing=1
fi

echo "=== building yelmox_bipolar ==="
make yelmox_bipolar || { echo "build failed"; exit 1; }

status=0
summary=()
result() {
    # result <PASS|FAIL> <test> <message>
    echo "  [$1] $2${3:+ ($3)}"
    summary+=("$1 $2")
    [ "$1" = "PASS" ] || status=1
}

log_has() {
    # log_has <log> <text>: the text in the log, with the lines joined (a long
    # list-directed write may wrap onto a new record with a leading blank).
    awk 'BEGIN { ORS = "" } { sub(/^ /, ""); print }' "$1" 2>/dev/null | grep -qF -- "$2"
}

run_test() {
    # run_test <test> [<expected error>]: <test> names the -p arrays, joined by
    # "+"; with an expected error the run must stop at init with it.
    local test="$1" expect="${2:-}" outdir="$OUT/${1//+/_}" name
    local -a params=() names=()
    IFS="+" read -r -a names <<< "$test"
    for name in "${names[@]}"; do
        if ! declare -p "$name" >/dev/null 2>&1; then
            result FAIL "$test" "no array $name"; return
        fi
        eval 'params+=("${'"$name"'[@]}")'
    done

    echo
    echo "=================================================================="
    echo "  $test${expect:+  (must stop: $expect)}"
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
        result FAIL "$test" "runme could not stage/launch the run"; return
    fi

    local log="$outdir/out.out" waited=0
    while true; do
        if grep -q "yelmox_bipolar: run complete" "$log" 2>/dev/null; then
            if [ -n "$expect" ]; then
                result FAIL "$test" "ran to completion, expected an error at init"
            else
                result PASS "$test"
            fi
            return
        fi
        if grep -qE "STOP [0-9]+|Fortran runtime error|Backtrace for|forrtl: severe" "$log" 2>/dev/null \
           || ! pgrep -f "yelmox_bipolar.x" >/dev/null 2>&1; then
            sleep 2                      # let the log flush
            if [ -n "$expect" ] && log_has "$log" "$expect"; then
                result PASS "$test" "stopped with the expected error"
            else
                result FAIL "$test" "stopped, see $log"
                grep -m3 -E ":: error|STOP [0-9]+|forrtl" "$log" | sed 's/^/         /'
            fi
            return
        fi
        sleep 10; waited=$((waited + 10))
        if [ "$waited" -ge "$TIMEOUT" ]; then
            result FAIL "$test" "timeout after ${TIMEOUT}s; still running, see $log"; return
        fi
    done
}

check_forcing() {
    # check_forcing <test> <run dir> [<reference run dir>]: the forcing of a
    # finished run (check_forcing.jl), against its continuous reference run.
    local test="$1"; shift
    if [ ! -f "$1/forcing_ts.nc" ]; then
        result FAIL "check $test" "no $1/forcing_ts.nc"; return
    fi
    if "$JULIA" --project=yelmox_bipolar yelmox_bipolar/check_forcing.jl "$@"; then
        result PASS "check $test"
    else
        result FAIL "check $test" "see the check output above"
    fi
}

if [ "$run_forcing" = "1" ]; then
    "$JULIA" --project=yelmox_bipolar -e 'using Pkg; Pkg.instantiate()' || exit 1
fi

for test in "${setups[@]}"; do
    run_test "$test"
    [ "$run_forcing" = "1" ] && check_forcing "$test" "$OUT/${test//+/_}"
done
if [ "$run_checks" = "1" ]; then
    for ((i = 0; i < ${#checks[@]}; i += 2)); do
        run_test "${checks[i]}" "${checks[i+1]}"
    done
fi
if [ "$run_forcing" = "1" ]; then
    # The last forcing test restarted from its bundle at 50 years.
    ref="$OUT/${test//+/_}"
    bip_restart50=("ctrl.restart=$(pwd)/$ref/restart-0.050-kyr" "ctrl.time_init=50")
    run_test "$test+bip_restart50"
    check_forcing "$test+bip_restart50" "$OUT/${test//+/_}_bip_restart50" "$ref"
fi

echo
echo "=== summary ==="
printf '  %s\n' "${summary[@]}"
exit "$status"
