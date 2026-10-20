#!/bin/bash
#
# Antarctic paleo transient, last glacial period (-130 kyr to +10 kyr, 32 km),
# restarting from the spin-up bundle at 15 kyr. Climate follows the glacial
# index input/alpha_combined_125kyr_interp.dat (index "a"). Base
# yelmox/yelmox_Antarctica.nml with antpal_physics and the forcing antpal_lgp
# of params.sh. Run from the yelmox root.

source "$(dirname "$0")/params.sh"

spinup_path=${1:-output/ant-paleo/spinup}
output_path=${2:-output/ant-paleo/lgp}

# Absolute path: the executable runs from inside the run dir, so a
# repo-root-relative path would not resolve.
restart="$(realpath -m "${spinup_path}")/restart-15.000-kyr"
if [ ! -d "${restart}" ]; then
    echo "run_lgp.sh: spin-up restart bundle not found: ${restart}" >&2
    exit 1
fi

run_params=(
    "ctrl.tstep_method=rel"     # forcing clock = model time (paleo sea level follows time)
    "ctrl.time_init=-130e3"     # [yr]
    "ctrl.time_end=10e3"        # [yr]
    "sim.relax=False"           # no relaxation
    "sim.opt=False"             # no optimization (cb_ref, tf_corr from the restart)
    "tm_1D.dt=100"              # [yr] 1D output
    "tm_2D.dt=1000"             # [yr] 2D output
    "tm_2Dsm.dt=1000000"        # [yr] small 2D output (off in practice)
    "tm_rst.dt=25e3"            # [yr] restart bundles
)

# Age tracing with elsa (Lagrangian isochrone layers), started fresh at
# time_init: the spin-up this run restarts from has no elsa layers.
elsa_params=(
    "ytrc.use_euler=False"                  # no Eulerian age tracer
    "ytrc.use_tracer=False"                 # no Lagrangian particles
    "ytrc.use_elsa=True"                    # Lagrangian layers (elsa)
    "ytrc.elsa_restart=False"               # start the layers fresh
    "ytrc.t_dep_source=elsa"                # deposition time from elsa
    "ytrc.time_end=10e3"                    # [yr] = ctrl.time_end (sizes elsa's layer stack)
    "ytrc.elsa_nml=input/forcing_db.nml"    # elsa group library
    "ytrc.elsa_group=elsa_ant_paleo"        # an isochrone every 1000 yr
)

# Submits to the queue by default; set runopts to change that (other queues, e.g. on
# albedo: runopts="-rs -q 48h -w 48:00:00", or -s to stage only).
runopts="${runopts:--rs -q shared -w 2-00:00:00 -m 10G}"

runme $runopts -e yelmox -n yelmox/yelmox_Antarctica.nml -o "${output_path}" \
      -p "${antpal_physics[@]}" "${antpal_lgp[@]}" "${run_params[@]}" "${elsa_params[@]}" \
         ctrl.restart="${restart}"
