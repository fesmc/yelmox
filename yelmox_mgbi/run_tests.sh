#!/usr/bin/env bash
#
# Smoke tests for the yelmox_mgbi bipolar driver (ocean-coupling port).
#
# Runs short bipolar simulations LOCALLY through the usual runme workflow: `-r`
# runs the executable directly (use `-rs` instead to submit to a queue). Each
# test writes to its own directory under tmp/. Run by hand with:
#
#     bash yelmox_mgbi/run_tests.sh              # tests 1 + 2
#     RUN_FWF=1 bash yelmox_mgbi/run_tests.sh    # also the freshwater-flux test
#
# Prerequisites:
#   - ice-sheet input data for GRL-16KM (north/Greenland) and ANT-32KM
#     (south/Antarctica) present under ice_data/.
#   - for Test 3 only: the hydrographic masks referenced by [ctrl]
#     hydro_mask_north/south (…_HYDRO-BASINS_zickfeld2004.nc).

cd "$(dirname "$0")/.." || exit 1        # repo root

NML="yelmox_mgbi/yelmox_mgbi_Bipolar.nml"
EXE="mgbi"
RUN_FWF="${RUN_FWF:-0}"

echo "=== building yelmox_mgbi ==="
make yelmox_mgbi || { echo "build failed"; exit 1; }

status=0
TIMEOUT="${TIMEOUT:-3600}"               # per-test wall-clock cap [s]
run_test() {
    # run_test <label> <outdir> <extra runme -p args...>
    local label="$1" outdir="$2"; shift 2
    echo
    echo "=================================================================="
    echo "  $label"
    echo "  -> $outdir"
    echo "=================================================================="
    rm -rf "$outdir"                     # fresh dir so the test is re-runnable

    # runme -r launches the executable in the BACKGROUND (…/out.out &) and returns
    # immediately, so poll the run log for a clean finish vs. an abort. Keep this
    # shell alive meanwhile, otherwise the backgrounded job is killed on exit.
    if ! runme -r -e "$EXE" -n "$NML" -o "$outdir" "$@"; then
        echo "  [FAIL] $label (runme could not stage/launch the run)"; status=1; return
    fi

    local log="$outdir/out.out" waited=0
    while true; do
        if grep -q "yelmox_mgbi: run complete" "$log" 2>/dev/null; then
            echo "  [PASS] $label (run complete)"; return
        fi
        if grep -qE "STOP [0-9]+|Fortran runtime error|Backtrace for" "$log" 2>/dev/null; then
            echo "  [FAIL] $label (aborted: $(grep -m1 -oE 'STOP [0-9]+|Fortran runtime error.*' "$log"))"
            echo "         Note: a Yelmo 'STOP 9' at time=0 is its own stability kill"
            echo "         (yelmo_check_kill), a property of the nml/data, not the driver."
            status=1; return
        fi
        if ! pgrep -f "yelmox_mgbi.x" >/dev/null 2>&1; then
            echo "  [FAIL] $label (process exited without a completion marker; see $log)"; status=1; return
        fi
        sleep 10; waited=$((waited + 10))
        if [ "$waited" -ge "$TIMEOUT" ]; then
            echo "  [FAIL] $label (timeout after ${TIMEOUT}s; still running, see $log)"; status=1; return
        fi
    done
}

# --- Test 1: baseline, ocean coupling OFF (two independent domains) -----------
# Confirms the refactor (shared bsl + explicit north/south) runs with no
# regression vs. the pre-port behaviour. Uses the nml defaults (active_obm=False).
run_test "Test 1: coupling OFF (baseline)" "tmp/mgbi_test1_baseline"

# --- Test 2: ocean box model ON; atmosphere <-> ocean <-> ice ----------------
# Exercises obm_update + atm2obm + obm2ism (box-model ocean temperature fed into
# marine_shelf). Freshwater flux (ism2obm) is OFF here because it needs the
# hydrographic mask files, which are not in every checkout (see Test 3).
run_test "Test 2: OBM on (atm2obm + obm2ism, no fwf)" "tmp/mgbi_test2_obm" \
    -p ctrl.active_obm=True ctrl.atm2obm=True ctrl.obm2ism=True

# --- Test 3: FULL coupling incl. freshwater flux (ism2obm) --------------------
# Needs the hydrographic masks referenced by [ctrl] hydro_mask_north/south.
# Disabled by default; run with RUN_FWF=1 once those files exist.
if [ "$RUN_FWF" = "1" ]; then
    run_test "Test 3: FULL coupling (+ ism2obm/fwf)" "tmp/mgbi_test3_fullcouple" \
        -p ctrl.active_obm=True ctrl.atm2obm=True ctrl.obm2ism=True \
           ctrl.ism2obm=True ctrl.couple_fwf_north=True ctrl.couple_fwf_south=True
else
    echo
    echo "=== Test 3 (full fwf coupling) skipped; set RUN_FWF=1 to enable (needs hydro masks) ==="
fi

echo
if [ "$status" = "0" ]; then
    echo "=== all requested mgbi smoke tests PASSED ==="
else
    echo "=== some mgbi smoke tests FAILED (see [FAIL] above) ==="
fi
exit "$status"
