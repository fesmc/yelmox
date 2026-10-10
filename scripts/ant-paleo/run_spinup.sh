#!/bin/bash
#
# Antarctic paleo spin-up (32 km): 15 kyr basal-friction and thermal-forcing
# optimization under present-day climate. Par file:
# yelmox/yelmox_Antarctica_paleo_spinup.nml. Run from the yelmox root.

output_path=${1:-output/ant-paleo/spinup}

# Submits to the queue by default; set runopts to change that (other queues, e.g. on
# albedo: runopts="-rs -q 48h -w 48:00:00", or -s to stage only).
runopts="${runopts:--rs -q shared -w 2-00:00:00 -m 10G}"

runme $runopts -e yelmox -n yelmox/yelmox_Antarctica_paleo_spinup.nml -o "${output_path}"
