#!/bin/bash
#
# Antarctic paleo spin-up (32 km): 15 kyr basal-friction and thermal-forcing
# optimization under present-day climate. Par file:
# yelmox/yelmox_Antarctica_paleo_spinup.nml. Run from the yelmox root.

output_path=${1:-output/ant-paleo/spinup}

runme -rs -q shared -e yelmox -w 2-00:00:00 -m 10G -n yelmox/yelmox_Antarctica_paleo_spinup.nml -o "${output_path}"
