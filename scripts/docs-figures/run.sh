#!/usr/bin/env bash
#
# Climate-only runs for the forcing figures of the docs pages (docs/figures/):
# each case runs its par file with the ice sheet off over a forcing time range,
# the climate updated and written every DT years (present-day topography, so
# the maps show the forcing only). plot.sh then draws the figures. See README.md.
#
# Usage (Levante, from anywhere; runs from the repo root):
#   scripts/docs-figures/run.sh [CASE ...]      (default: all cases)
#   scripts/docs-figures/run.sh --list          the cases as presets of the forcing
#       builder (scripts/forcing-builder/export.jl): one line per case, tab
#       separated, case, par file, then the -p overrides of the experiment and
#       its timeline (not the climate-only ones)
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

# Shared -p arrays of the experiments, referenced as @<array> in the
# overrides column.
source scripts/grl-paleo/params.sh
source scripts/ant-paleo/params.sh
source scripts/ismip7/params.sh
source scripts/1pctCO2/params.sh
source scripts/tipmip/params.sh

# case | par file | timeline group | tstep_method | t0 | t1 | dt | extra overrides
# (@name: the elements of the bash array name)
# (rel: time relative to 1950; cal: calendar years; const: the forcing at t0)
CASES=$(cat <<'EOF'
grl_snap   |yelmox/yelmox_Greenland.nml               |ctrl      |rel  |-130e3 |0      |1000 |
grl_ramp   |yelmox/yelmox_Greenland.nml               |ctrl      |rel  |0      |2000   |500  |domain.grid_hub=GRL-32KM domain.grid_isos=GRL-32KM domain.grid_clim=GRL-32KM domain.grid_surface=GRL-32KM domain.grid_shelf=GRL-32KM regions.basin_sets=Zwally2012 forcing.series==orbital,dTo,ramp forcing.terms==orbital,ocean,ramp
ant_pd     |yelmox/yelmox_Antarctica.nml              |ctrl      |const|0      |10     |10   |
ant_paleo  |yelmox/yelmox_Antarctica.nml              |ctrl      |rel  |-130e3 |0      |1000 |@antpal_physics @antpal_lgp sim.relax=False sim.opt=False ytill.method=1
grlpal_pd  |yelmox/yelmox_GreenlandPaleo.nml          |ctrl      |const|0      |10     |10   |
grlpal_a   |yelmox/yelmox_GreenlandPaleo.nml          |ctrl      |rel  |-130e3 |0      |1000 |@grlpal_lgp
grlpal_b   |yelmox/yelmox_GreenlandPaleo.nml          |ctrl      |rel  |-21e3  |0      |500  |@grlpal_recon
lis_lgm    |yelmox/yelmox_LIS.nml                     |ctrl      |const|-21e3  |10     |10   |
lis_pd     |yelmox/yelmox_LIS.nml                     |ctrl      |const|0      |10     |10   |series_glacial.f_min=0
north_lgm  |yelmox/yelmox_North.nml                   |ctrl      |const|-21e3  |10     |10   |
north_pd   |yelmox/yelmox_North.nml                   |ctrl      |const|0      |10     |10   |series_glacial.f_min=0
pyrenees   |yelmox/yelmox_Pyrenees.nml                |ctrl      |rel  |-21e3  |0      |1000 |
srg_pd     |yelmox/yelmox_SRG.nml                     |ctrl      |const|1950   |10     |10   |
ismip7_grl |yelmox/yelmox_esm_Greenland.nml           |transient |cal  |1975   |2300   |25   |ctrl.run_step=transient esm.experiment=ssp585 esm.use_esm=True esm.use_hist=True esm.use_proj=True
ismip7_ant |yelmox/yelmox_esm_Antarctica.nml          |transient |cal  |1975   |2300   |25   |@ismip7_ant ctrl.run_step=transient esm.experiment=ssp585 esm.use_esm=True esm.use_hist=True esm.use_proj=True
pct_grl    |yelmox/yelmox_esm_Greenland.nml           |transient |cal  |2020   |2160   |20   |@pct_grl ctrl.run_step=transient esm.experiment=1pctCO2 esm.use_esm=True esm.use_hist=False esm.use_proj=True
pct_ant    |yelmox/yelmox_esm_Antarctica.nml          |transient |cal  |2020   |2160   |20   |@ismip7_ant @pct_ant ctrl.run_step=transient esm.experiment=1pctCO2 esm.use_esm=True esm.use_hist=False esm.use_proj=True
tipmip_grl |yelmox/yelmox_esm_Greenland.nml           |transient |cal  |0      |230    |10   |@tipmip_grl ctrl.run_step=transient esm.experiment=esm-up2p0 esm.use_esm=True esm.use_hist=False esm.use_proj=True
EOF
)

LIST=0; [ "${1:-}" = "--list" ] && { LIST=1; shift; }

want () { [ $# -eq 0 ] && return 0; local c; for c in "${SEL[@]}"; do [ "$c" = "$1" ] && return 0; done; return 1; }
SEL=("$@")

while IFS='|' read -r name par g ts t0 t1 dt extra; do
    name=$(echo $name); par=$(echo $par); g=$(echo $g); ts=$(echo $ts)
    t0=$(echo $t0); t1=$(echo $t1); dt=$(echo $dt)
    [ ${#SEL[@]} -gt 0 ] && ! want "$name" && continue
    # The timeline of the case.
    case $ts in
        const) tl=("$g.tstep_method=const" "$g.tstep_const=$t0" "$g.time_init=0" "$g.time_end=$t1") ;;
        rel)   tl=("$g.tstep_method=rel" "$g.time_init=$t0" "$g.time_end=$t1") ;;
        cal)   tl=("$g.tstep_method=cal" "$g.time_init=$t0" "$g.time_end=$t1") ;;
    esac
    read -r -a toks <<< "$extra"
    ex=()
    for tok in ${toks[@]+"${toks[@]}"}; do
        if [[ $tok == @* ]]; then ref="${tok#@}[@]"; ex+=("${!ref}"); else ex+=("$tok"); fi
    done
    if [ $LIST = 1 ]; then
        (IFS=$'\t'; echo "${name}"$'\t'"${par}"$'\t'"${ex[*]+${ex[*]}}${ex[*]+$'\t'}${tl[*]}")
        continue
    fi
    # Climate-only: ice sheet and isostasy off, climate updated and written
    # every dt (write_clim on: the ESM par files have it off), the final
    # restart bundle only.
    o=("comps.with_ice_sheet=False" "comps.with_isostasy=False" "comps.dt_clim=$dt" "$g.dtt=$dt"
       "tm_2D.method=const" "tm_2D.dt=$dt" "tm_2Dsm.method=const" "tm_2Dsm.dt=$dt"
       "tm_1D.method=const" "tm_1D.dt=$dt" "tm_rst.method=none" "output.write_clim=True" "${tl[@]}")
    echo ">> $name: $par -> $OUT/$name"
    # The case's overrides first, the climate-only ones last: they win.
    runme $MODE "${RUNOPTS_ARR[@]}" --jobname "docfig-$name" -e yelmox -n "$par" \
        -o "$OUT/$name" -p ${ex[@]+"${ex[@]}"} "${o[@]}" > /dev/null
done <<< "$CASES"
