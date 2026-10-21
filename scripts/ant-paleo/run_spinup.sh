#!/bin/bash
#
# Antarctic paleo spin-up (32 km): 15 kyr basal-friction and thermal-forcing
# optimization under present-day climate. Base yelmox/yelmox_Antarctica.nml
# with the 32 km setup antpal_physics of params.sh. Run from the yelmox root.

source "$(dirname "$0")/params.sh"

output_path=${1:-output/ant-paleo/spinup}

spinup_params=(
    "series_dTo.factor=0.40"    # [1] ocean anomaly: fraction of the air temperature anomaly
    "tm_1D.dt=100"              # [yr] 1D output
    "tm_2D.dt=5000"             # [yr] 2D output
    "tm_2Dsm.dt=300e3"          # [yr] small 2D output (off in practice)
    "tm_rst.dt=15e3"            # [yr] restart bundles (the transient restarts from 15 kyr)
)

# Submits to the queue by default; set runopts to change that (other queues, e.g. on
# albedo: runopts="-rs -q 48h -w 48:00:00", or -s to stage only).
runopts="${runopts:--rs -q shared -w 2-00:00:00 -m 10G}"

runme $runopts -e yelmox -n yelmox/yelmox_Antarctica.nml -o "${output_path}" \
      -p "${antpal_physics[@]}" "${spinup_params[@]}"
