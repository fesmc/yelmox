#!/bin/bash
#
# Greenland paleo spin-up (GRL-PAL, 16 km): 15 kyr under present-day climate,
# relaxation early on, no optimization. Base yelmox/yelmox_GreenlandPaleo.nml.
# Run from the yelmox root.

output_path=${1:-output/grl-paleo/spinup}

spinup_params=(
    "sim.relax=True"        # relax the topography towards the reference ([relax])
    "tm_2D.dt=1000"         # [yr] 2D output
    "tm_rst.dt=15e3"        # [yr] restart bundles (the transient restarts from 15 kyr)
)

# Submits to the queue by default; set runopts to change that (other queues, e.g. on
# albedo: runopts="-rs -q 48h -w 48:00:00", or -s to stage only).
runopts="${runopts:--rs -q shared -w 2-00:00:00 -m 10G}"

runme $runopts -e yelmox -n yelmox/yelmox_GreenlandPaleo.nml -o "${output_path}" \
      -p "${spinup_params[@]}"
