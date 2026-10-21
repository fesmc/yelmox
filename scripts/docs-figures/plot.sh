#!/usr/bin/env bash
#
# Forcing figures of the docs pages (docs/figures/<case>.png) from the runs of
# run.sh. Usage (from anywhere): scripts/docs-figures/plot.sh [CASE ...]
#
# Environment: OUT (the output root of run.sh, default output/docs-figures),
# FIG (the figure folder, default docs/figures).
#
set -euo pipefail
cd "$(dirname "$0")/../.."

OUT=${OUT:-output/docs-figures}
FIG=${FIG:-docs/figures}
mkdir -p "$FIG"
J="julia --project=scripts/docs-figures scripts/docs-figures/plot.jl"
SEL=("$@")
want () { [ ${#SEL[@]} -eq 0 ] && return 0; local c; for c in "${SEL[@]}"; do [ "$c" = "$1" ] && return 0; done; return 1; }

# figure | run | plot.jl options
while IFS='|' read -r fig run opts; do
    fig=$(echo $fig); run=$(echo $run)
    want "$fig" || continue
    eval "$J $OUT/$run $FIG/$fig.png $opts"
done <<'EOF2'
grl_snap   | grl_snap   | --ref 0 --at -21000 -125000 --title "yelmox_Greenland.nml: snapshots on the glacial index"
grl_ramp   | grl_ramp   | --ref 0 --at 2000 --title "yelmox_Greenland.nml + term_ramp: air temperature ramp"
ant_pd     | ant_pd     | --ref 10 --abs --title "yelmox_Antarctica.nml: present-day reference climate"
ant_paleo  | ant_paleo  | --ref 0 --at -21000 -125000 --title "yelmox_Antarctica.nml + antpal_physics, antpal_lgp"
grlpal_pd  | grlpal_pd  | --ref 10 --abs --title "yelmox_GreenlandPaleo.nml: ERA5 + MAR - ERA5 bias"
grlpal_a   | grlpal_a   | --ref 0 --at -21000 -125000 --title "yelmox_GreenlandPaleo.nml + grlpal_lgp (A)"
grlpal_b   | grlpal_b   | --ref 0 --at -21000 -12500 --title "yelmox_GreenlandPaleo.nml + grlpal_recon (B)"
lis        | lis_lgm    | --ref 10 --at 10 --ref-run $OUT/lis_pd --title "yelmox_LIS.nml: LGM (index a = 1) minus present day (a = 0)"
north      | north_lgm  | --ref 10 --at 10 --ref-run $OUT/north_pd --title "yelmox_North.nml: LGM (index a = 1) minus present day (a = 0)"
pyrenees   | pyrenees   | --ref 0 --at -21000 --title "yelmox_Pyrenees.nml"
srg_pd     | srg_pd     | --ref 10 --abs --title "yelmox_SRG.nml: present-day reference climate"
ismip7_grl | ismip7_grl | --ref 1975 --at 2100 2300 --title "yelmox_esm_Greenland.nml: ssp585"
ismip7_ant | ismip7_ant | --ref 1975 --at 2100 2300 --title "yelmox_esm_Antarctica.nml + ismip7_ant: ssp585"
pct_grl    | pct_grl    | --ref 2020 --at 2160 --title "yelmox_esm_Greenland.nml + pct_grl"
pct_ant    | pct_ant    | --ref 2020 --at 2160 --title "yelmox_esm_Antarctica.nml + ismip7_ant, pct_ant"
tipmip_grl | tipmip_grl | --ref 0 --at 230 --title "yelmox_esm_Greenland.nml + tipmip_grl: esm-up2p0"
EOF2
