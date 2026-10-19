#!/bin/bash
#
# Greenland LGM spin-up (GRL-PAL, 16 km) for the deglaciation: cold start from
# GLAC-1D at 21 ka, 20 kyr with the forcing held at 21 ka. Par file:
# yelmox/yelmox_GreenlandPaleo_lgm.nml. Run from the yelmox root.

output_path=${1:-output/grl-paleo/lgm}

# Submits to the queue by default; set runopts to change that (other queues, e.g. on
# albedo: runopts="-rs -q 48h -w 48:00:00", or -s to stage only).
runopts="${runopts:--rs -q shared -w 2-00:00:00 -m 10G}"

runme $runopts -e yelmox -n yelmox/yelmox_GreenlandPaleo_lgm.nml -o "${output_path}"
