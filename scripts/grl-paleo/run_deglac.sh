#!/bin/bash
#
# Greenland deglaciation after Tabone et al. (2024), -21 ka to 0 (GRL-PAL,
# 16 km), restarting from the LGM spin-up bundle at 20 kyr. Par file:
# yelmox/yelmox_GreenlandPaleo_deglac.nml. Run from the yelmox root.

lgm_path=${1:-output/grl-paleo/lgm}
output_path=${2:-output/grl-paleo/deglac}

# Absolute path: the executable runs from inside the run dir, so a
# repo-root-relative path would not resolve.
restart="$(realpath -m "${lgm_path}")/restart-20.000-kyr"
if [ ! -d "${restart}" ]; then
    echo "run_deglac.sh: LGM spin-up restart bundle not found: ${restart}" >&2
    exit 1
fi

# Submits to the queue by default; set runopts to change that (other queues, e.g. on
# albedo: runopts="-rs -q 48h -w 48:00:00", or -s to stage only).
runopts="${runopts:--rs -q shared -w 2-00:00:00 -m 10G}"

runme $runopts -e yelmox -n yelmox/yelmox_GreenlandPaleo_deglac.nml -o "${output_path}" \
      -p ctrl.restart="${restart}"
