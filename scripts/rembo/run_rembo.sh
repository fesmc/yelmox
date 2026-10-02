#!/usr/bin/env bash
#
# Default REMBO-coupled Greenland run (yelmox with climate = "rembo").
#
# Single grid (GRL-16KM): ice sheet + isostasy + REMBO atmosphere/SMB + snapclim
# ocean + marine_shelf, present-day style (no hysteresis, no optimization). REMBO
# provides the surface mass balance; the shared multigrid couplers land it (and
# isostasy / marine_shelf) on the Yelmo grid each step. See
# yelmox/yelmox_rembo_Greenland.nml for the full configuration, and
# yelmox/rembo_Greenland.nml for REMBO's own parameters (staged into the
# run directory automatically via .runme/info.json).
#
# Usage (from anywhere):
#     scripts/rembo/run_rembo.sh [output_dir]
#
# Runs locally by default. To submit to the queue instead:
#     runopts='-rs -q 12h -w 10:00:00' scripts/rembo/run_rembo.sh
#
# Build first with:  make yelmox rembo=1

cd "$(dirname "$0")/../.." || exit 1     # repo root

EXE="yelmox"
NML="yelmox/yelmox_rembo_Greenland.nml"
OUT="${1:-output/rembo}"

runopts="${runopts:--r}"                 # local run by default; set runopts='-rs ...' to submit

runme $runopts -e "$EXE" -n "$NML" -o "$OUT"
