#!/usr/bin/env bash
#
# Antarctica multigrid resolution test -- forward projections (yelmox).
#
# Follow-up to run_mg_resolution.sh. Branches 1 kyr forward runs off the
# spun-up (t=25 kyr) state of each COMPLETED core/multigrid run, under two
# forcings:
#   - ctrl  : present-day climate held fixed (dTa=dTo=0);
#   - warm  : uniform +2 K atmosphere and +2 K ocean anomaly (dTa=dTo=+2).
#
# Only the three Yelmo-32KM spin-ups finished within the wall-clock limit, so
# the projections branch from those (marine_shelf grid = 8 / 16 / 32 KM). This
# isolates how the marine_shelf resolution changes the warming response -- the
# ocean anomaly drives shelf melt (bmb_shlf), the field the finer mshlf grid
# resolves, so the divergence (if any) should appear here rather than in the
# present-day spin-up.
#
# Each run restarts from the parent's restart bundle folder via
# coupling.restart, with:
#   - coupling.equil_method=none : freeze the optimized cf_ref / tf_corr that
#       were loaded from the restart (no further optimization);
#   - yelmo.restart_relax=0      : free evolution from the spun-up state (do NOT
#       relax H_ice back toward the input PD topography). *** Sanity-check this:
#       set >0 only if a short relaxation is wanted for numerical stability. ***
# The grids for each projection match its parent spin-up exactly.
#
# Usage (from anywhere):
#     yelmox/run_mg_forward.sh
#
# Submits to the queue by default. Set runopts='-r' to run locally instead.

ROOT="$(cd "$(dirname "$0")/.." && pwd)" || exit 1   # repo root (absolute)
cd "$ROOT" || exit 1

runopts='-rs -q compute -w 04:00:00'

EXE="yelmox"
NML="yelmox/yelmox_Antarctica.nml"
SPINROOT="$ROOT/output/mg_opt"     # where the completed spin-ups live
OUTROOT="output/mg_fwd"            # where projections are written
RESTART_KYR="restart-25.000-kyr"   # spun-up bundle to branch from
TIME_END=1000                      # [yr] projection length

# run_fwd <parent_subdir> <yelmo_grid> <mshlf_grid> <scenario> <dTa> <dTo>
run_fwd() {
    local parent="$1" ygrid="$2" mgrid="$3" scen="$4" dTa="$5" dTo="$6"
    local restart="$SPINROOT/$parent/$RESTART_KYR"
    local out="$OUTROOT/${parent}_${scen}"

    if [ ! -f "$restart/yelmo_restart.nc" ]; then
        echo "SKIP $out : missing restart bundle $restart" >&2
        return
    fi

    runme $runopts -e "$EXE" -n "$NML" -o "$out" \
        -p yelmo.grid_name="$ygrid" htopo.grid_name="$mgrid" \
           coupling.grid_mshlf="$mgrid" coupling.grid_isos="$ygrid" \
           coupling.grid_clim=ANT-32KM coupling.grid_smb="$ygrid" \
           coupling.restart="$restart" coupling.equil_method=none \
           yelmo.restart_relax=0 \
           ctrl.time_init=0 ctrl.time_end="$TIME_END" \
           tm_2D.dt=100 ytopo.gz_Hg1=100 \
           snap.dTa_const="$dTa" snap.dTo_const="$dTo"
}

# For each completed spin-up: a ctrl run and a +2/+2 warming run.
# run_fwd <parent>      <yelmo>   <mshlf>   <scen> <dTa> <dTo>
run_fwd y32KM_m32KM  ANT-32KM  ANT-32KM  ctrl   0  0
run_fwd y32KM_m32KM  ANT-32KM  ANT-32KM  warm   2  2
run_fwd y32KM_m16KM  ANT-32KM  ANT-16KM  ctrl   0  0
run_fwd y32KM_m16KM  ANT-32KM  ANT-16KM  warm   2  2
run_fwd y32KM_m8KM   ANT-32KM  ANT-8KM   ctrl   0  0
run_fwd y32KM_m8KM   ANT-32KM  ANT-8KM   warm   2  2
