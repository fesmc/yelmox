#!/bin/bash
#
# Greenland paleo spin-up (GRL-PAL, 16 km): 15 kyr under present-day climate,
# no optimization. Par file:
# yelmox/yelmox_GreenlandPaleo_spinup.nml. Run from the yelmox root.

output_path=${1:-output/grl-paleo/spinup}

# Submits to the queue by default; set runopts to change that (other queues, e.g. on
# albedo: runopts="-rs -q 48h -w 48:00:00", or -s to stage only).
runopts="${runopts:--rs -q shared -w 2-00:00:00 -m 10G}"

runme $runopts -e yelmox -n yelmox/yelmox_GreenlandPaleo_spinup.nml -o "${output_path}"
