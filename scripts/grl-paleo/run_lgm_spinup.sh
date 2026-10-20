#!/bin/bash
#
# Greenland LGM spin-up (GRL-PAL, 16 km) for the deglaciation: cold start from
# GLAC-1D at 21 ka, 20 kyr with the forcing held at 21 ka (constant time -21e3:
# reconstructions, insolation, sea level). Forcing grlpal_recon of params.sh.
# Base yelmox/yelmox_GreenlandPaleo.nml. Run from the yelmox root.

source "$(dirname "$0")/params.sh"

output_path=${1:-output/grl-paleo/lgm}

run_params=(
    "ctrl.tstep_const=-21e3"    # [yr] forcing held at 21 ka
    "ctrl.time_end=20e3"        # [yr]
    "tm_1D.dt=10"               # [yr] 1D output
    "tm_2D.dt=1000"             # [yr] 2D output
)

# Cold start from the GLAC-1D ice at 21 ka.
init_params=(
    "sim.init_method=recon"     # ice from recon_path where recon_regions
    "sim.recon_path=fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-GLAC1D.nc"
    "sim.recon_var=H_ice"       # ice thickness in recon_path
    "sim.recon_time=-21e3"      # [yr] nearest slice of recon_path
    "sim.recon_regions=all"     # imposed everywhere
)

# Submits to the queue by default; set runopts to change that (other queues, e.g. on
# albedo: runopts="-rs -q 48h -w 48:00:00", or -s to stage only).
runopts="${runopts:--rs -q shared -w 2-00:00:00 -m 10G}"

runme $runopts -e yelmox -n yelmox/yelmox_GreenlandPaleo.nml -o "${output_path}" \
      -p "${grlpal_recon[@]}" "${run_params[@]}" "${init_params[@]}"
