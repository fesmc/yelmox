#!/usr/bin/env bash
#
# Climate-only runs for the forcing figures of the docs pages (docs/figures/):
# each case runs its par file with the ice sheet off over a forcing time range,
# the climate updated and written every DT years (present-day topography, so
# the maps show the forcing only). plot.sh then draws the figures. See README.md.
#
# Usage (Levante, from anywhere; runs from the repo root):
#   scripts/docs-figures/run.sh [CASE ...]      (default: all cases)
#
# Environment:
#   OUT       output root (default output/docs-figures)
#   RUNOPTS   slurm options for runme (default "-q compute -w 00:30:00 --omp 16")
#   STAGE=1   stage the runs only (submit scripts written, not submitted)
#
set -euo pipefail
cd "$(dirname "$0")/../.."     # repo root: runme needs .runme/ and the data links

OUT=${OUT:-output/docs-figures}
RUNOPTS="${RUNOPTS:--q compute -w 00:30:00 --omp 16}"
read -r -a RUNOPTS_ARR <<< "$RUNOPTS"
MODE="-rs"; [ "${STAGE:-0}" = 1 ] && MODE="-s"

# case | par file | timeline group | tstep_method | t0 | t1 | dt | extra overrides
# (rel: time relative to 1950; cal: calendar years; const: the forcing at t0)
CASES=$(cat <<'EOF'
grl_snap   |yelmox/yelmox_Greenland.nml               |ctrl      |rel  |-130e3 |0      |1000 |
grl_ramp   |yelmox/yelmox_Greenland_ramp.nml          |ctrl      |rel  |0      |2000   |500  |
ant_pd     |yelmox/yelmox_Antarctica.nml              |ctrl      |const|0      |10     |10   |
ant_paleo  |yelmox/yelmox_Antarctica_paleo_lgp.nml    |ctrl      |rel  |-130e3 |0      |1000 |
grlpal_pd  |yelmox/yelmox_GreenlandPaleo.nml          |ctrl      |const|0      |10     |10   |
grlpal_a   |yelmox/yelmox_GreenlandPaleo_lgp.nml      |ctrl      |rel  |-130e3 |0      |1000 |
grlpal_b   |yelmox/yelmox_GreenlandPaleo_deglac.nml   |ctrl      |rel  |-21e3  |0      |500  |
lis_lgm    |yelmox/yelmox_LIS.nml                     |ctrl      |const|-21e3  |10     |10   |
lis_pd     |yelmox/yelmox_LIS.nml                     |ctrl      |const|0      |10     |10   |idx_a.f_min=0
north_lgm  |yelmox/yelmox_North.nml                   |ctrl      |const|-21e3  |10     |10   |
north_pd   |yelmox/yelmox_North.nml                   |ctrl      |const|0      |10     |10   |idx_a.f_min=0
pyrenees   |yelmox/yelmox_Pyrenees.nml                |ctrl      |rel  |-120e3 |0      |1000 |
srg_pd     |yelmox/yelmox_SRG.nml                     |ctrl      |const|1950   |10     |10   |
ismip7_grl |yelmox/yelmox_esm_Greenland.nml           |transient |cal  |1975   |2300   |25   |ctrl.run_step=transient esm.experiment=ssp585 esm.use_esm=True esm.use_hist=True esm.use_proj=True
ismip7_ant |yelmox/yelmox_esm_Antarctica_ismip7.nml   |transient |cal  |1975   |2300   |25   |ctrl.run_step=transient esm.experiment=ssp585 esm.use_esm=True esm.use_hist=True esm.use_proj=True
pct_grl    |yelmox/yelmox_esm_Greenland_1pctCO2.nml   |transient |cal  |2020   |2160   |20   |ctrl.run_step=transient esm.experiment=1pctCO2 esm.use_esm=True esm.use_hist=False esm.use_proj=True
pct_ant    |yelmox/yelmox_esm_Antarctica_1pctCO2.nml  |transient |cal  |2020   |2160   |20   |ctrl.run_step=transient esm.experiment=1pctCO2 esm.use_esm=True esm.use_hist=False esm.use_proj=True
tipmip_grl |yelmox/yelmox_esm_Greenland_tipmip.nml    |transient |cal  |0      |230    |10   |ctrl.run_step=transient esm.experiment=esm-up2p0 esm.use_esm=True esm.use_hist=False esm.use_proj=True
EOF
)

want () { [ $# -eq 0 ] && return 0; local c; for c in "${SEL[@]}"; do [ "$c" = "$1" ] && return 0; done; return 1; }
SEL=("$@")

while IFS='|' read -r name par g ts t0 t1 dt extra; do
    name=$(echo $name); par=$(echo $par); g=$(echo $g); ts=$(echo $ts)
    t0=$(echo $t0); t1=$(echo $t1); dt=$(echo $dt)
    [ ${#SEL[@]} -gt 0 ] && ! want "$name" && continue
    # Climate-only: ice sheet and isostasy off, climate updated and written
    # every dt, the final restart bundle only.
    o=("comps.with_ice_sheet=False" "comps.with_isostasy=False" "comps.dt_clim=$dt" "$g.dtt=$dt"
       "tm_2D.method=const" "tm_2D.dt=$dt" "tm_2Dsm.method=const" "tm_2Dsm.dt=$dt"
       "tm_1D.method=const" "tm_1D.dt=$dt" "tm_rst.method=none")
    case $ts in
        const) o+=("$g.tstep_method=const" "$g.tstep_const=$t0" "$g.time_init=0" "$g.time_end=$t1") ;;
        rel)   o+=("$g.tstep_method=rel" "$g.time_init=$t0" "$g.time_end=$t1") ;;
        cal)   o+=("$g.tstep_method=cal" "$g.time_init=$t0" "$g.time_end=$t1") ;;
    esac
    read -r -a ex <<< "$extra"
    echo ">> $name: $par -> $OUT/$name"
    runme $MODE "${RUNOPTS_ARR[@]}" --jobname "docfig-$name" -e yelmox -n "$par" \
        -o "$OUT/$name" -p "${o[@]}" ${ex[@]+"${ex[@]}"} > /dev/null
done <<< "$CASES"
