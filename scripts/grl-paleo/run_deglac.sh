#!/bin/bash
#
# Greenland deglaciation after Tabone et al. (2024), -21 ka to 0 (GRL-PAL,
# 16 km), restarting from the LGM spin-up bundle at 20 kyr (run_spinup.sh recon
# -21e3): the reconstructions in time (forcing grlpal_recon of params.sh),
# transient insolation, sea level from file, elsa isochrones. Base
# yelmox/yelmox_GreenlandPaleo.nml. Run from the yelmox root.

source "$(dirname "$0")/params.sh"

lgm_path=${1:-output/grl-paleo/spinup-recon-21ka}
output_path=${2:-output/grl-paleo/deglac}

# Absolute path: the executable runs from inside the run dir, so a
# repo-root-relative path would not resolve.
restart="$(realpath -m "${lgm_path}")/restart-20.000-kyr"
if [ ! -d "${restart}" ]; then
    echo "run_deglac.sh: LGM spin-up restart bundle not found: ${restart}" >&2
    exit 1
fi

run_params=(
    "ctrl.tstep_method=rel"     # forcing time = model time
    "ctrl.time_init=-21e3"      # [yr]
    "ctrl.time_end=0.0"         # [yr]
    "comps.dt_clim=100"         # [yr] the reconstructions are read and remapped at each update
    "tm_1D.dt=10"               # [yr] 1D output
    "tm_2D.dt=500"              # [yr] 2D output
    "tm_2Dsm.dt=100"            # [yr] small 2D output
    "tm_rst.dt=7e3"             # [yr] restart bundles
)

# Age tracing with elsa (Lagrangian isochrone layers), started fresh at
# time_init: the LGM spin-up this run restarts from has no elsa layers.
elsa_params=(
    "ytrc.use_euler=False"                  # no Eulerian age tracer
    "ytrc.use_tracer=False"                 # no Lagrangian particles
    "ytrc.use_elsa=True"                    # Lagrangian layers (elsa)
    "ytrc.elsa_restart=False"               # start the layers fresh
    "ytrc.t_dep_source=elsa"                # deposition time from elsa
    "ytrc.time_end=0.0"                     # [yr] = ctrl.time_end (sizes elsa's layer stack)
    "ytrc.elsa_nml=input/forcing_db.nml"    # elsa group library
    "ytrc.elsa_group=elsa_grl_deglac"       # an isochrone every 500 yr
)

# Submits to the queue by default; set runopts to change that (other queues, e.g. on
# albedo: runopts="-rs -q 48h -w 48:00:00", or -s to stage only).
runopts="${runopts:--rs -q shared -w 2-00:00:00 -m 10G}"

runme $runopts -e yelmox -n yelmox/yelmox_GreenlandPaleo.nml -o "${output_path}" \
      -p "${grlpal_recon[@]}" "${run_params[@]}" "${elsa_params[@]}" \
         ctrl.restart="${restart}"
