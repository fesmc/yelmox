#!/bin/bash
#
# Greenland paleo spin-up (GRL-PAL, 16 km): 20 kyr with the forcing held at a
# constant time, no optimization. Base yelmox/yelmox_GreenlandPaleo.nml. Run
# from the yelmox root.
#
#   run_spinup.sh <forcing> <time> [output_path] [extra -p values ...]
#
#   forcing  pd (reference climate only), lgp (A) or recon (B): the array
#            grlpal_<forcing> of params.sh
#   time     [yr] forcing time (ctrl.tstep_const): climate, insolation, sea
#            level; e.g. 0 (present day) or -21e3 (LGM)
#
# The ice starts from the present-day topography at time 0 and from GLAC-1D at
# <time> otherwise. Output goes to output/grl-paleo/spinup-<forcing>-<ka>ka by
# default. Examples:
#   run_spinup.sh pd 0 output/grl-paleo/spinup sim.relax=True    (A)
#   run_spinup.sh recon -21e3 output/grl-paleo/lgm               (B)

source "$(dirname "$0")/params.sh"

if [ $# -lt 2 ]; then
    echo "usage: run_spinup.sh <forcing: pd|lgp|recon> <time [yr]> [output_path] [extra -p ...]" >&2
    exit 1
fi
forcing=$1
time=$2
shift 2

if ! declare -p "grlpal_${forcing}" > /dev/null 2>&1; then
    echo "run_spinup.sh: no forcing array grlpal_${forcing} in params.sh" >&2
    exit 1
fi
ref="grlpal_${forcing}[@]"
forcing_params=("${!ref}")

ka=$(awk -v t="${time}" 'BEGIN { printf "%g", -t/1000 }')
output_path=${1:-output/grl-paleo/spinup-${forcing}-${ka}ka}
[ $# -gt 0 ] && shift

run_params=(
    "ctrl.tstep_const=${time}"  # [yr] forcing held at this time
    "ctrl.time_end=20e3"        # [yr]
    "tm_1D.dt=10"               # [yr] 1D output
    "tm_2D.dt=1000"             # [yr] 2D output
    "tm_rst.dt=20e3"            # [yr] restart bundles (the transients restart from 20 kyr)
)

# Cold start from the GLAC-1D ice at <time> (the present-day topography at 0).
init_params=()
if [ "${ka}" != "0" ]; then
    init_params=(
        "sim.init_method=recon"     # ice from recon_path where recon_regions
        "sim.recon_path=fesmdata/{domain}/{grid_name}/{grid_name}_IceHistory-GLAC1D.nc"
        "sim.recon_var=H_ice"       # ice thickness in recon_path
        "sim.recon_time=${time}"    # [yr] nearest slice of recon_path
        "sim.recon_regions=all"     # imposed everywhere
    )
fi

# Submits to the queue by default; set runopts to change that (other queues, e.g. on
# albedo: runopts="-rs -q 48h -w 48:00:00", or -s to stage only).
runopts="${runopts:--rs -q shared -w 2-00:00:00 -m 10G}"

runme $runopts -e yelmox -n yelmox/yelmox_GreenlandPaleo.nml -o "${output_path}" \
      -p "${forcing_params[@]}" "${run_params[@]}" "${init_params[@]}" "$@"
