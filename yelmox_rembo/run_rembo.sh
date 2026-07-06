#!/usr/bin/env bash
#
# Default REMBO-coupled Greenland run (yelmox_rembo).
#
# Single grid (GRL-16KM): ice sheet + isostasy + REMBO atmosphere/SMB + snapclim
# ocean + marine_shelf, present-day style (no hysteresis, no optimization). REMBO
# provides the surface mass balance; the shared multigrid couplers land it (and
# isostasy / marine_shelf) on the Yelmo grid each step. See
# yelmox_rembo/yelmox_rembo_Greenland.nml for the full configuration, and
# yelmox_rembo/rembo_Greenland.nml for REMBO's own parameters (staged into the
# run directory automatically via .runme/info.json).
#
# Usage (from anywhere):
#     yelmox_rembo/run_rembo.sh [output_dir]
#
# Runs locally by default. To submit to the queue instead:
#     runopts='-rs -q 12h -w 10:00:00' yelmox_rembo/run_rembo.sh
#
# Build first with:  make yelmox_rembo

cd "$(dirname "$0")/.." || exit 1        # repo root

EXE="rembo"
NML="yelmox_rembo/yelmox_rembo_Greenland.nml"
OUT="${1:-output/rembo}"

runopts="${runopts:--r}"                 # local run by default; set runopts='-rs ...' to submit

runme $runopts -e "$EXE" -n "$NML" -o "$OUT"
