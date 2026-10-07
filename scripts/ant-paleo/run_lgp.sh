#!/bin/bash
#
# Antarctic paleo transient, last glacial period (-130 kyr to +10 kyr, 32 km),
# restarting from the spin-up bundle at 15 kyr. Climate follows the glacial
# index input/alpha_combined_125kyr_interp.dat (snapclim "snap_1ind"). Par file:
# yelmox/yelmox_Antarctica_paleo_lgp.nml. Run from the yelmox root.

spinup_path=${1:-output/ant-paleo/spinup}
output_path=${2:-output/ant-paleo/lgp}

# Absolute path: the executable runs from inside the run dir, so a
# repo-root-relative path would not resolve.
restart="$(realpath -m "${spinup_path}")/restart-15.000-kyr"
if [ ! -d "${restart}" ]; then
    echo "run_lgp.sh: spin-up restart bundle not found: ${restart}" >&2
    exit 1
fi

runme -rs -q shared -e yelmox -w 2-00:00:00 -m 10G -n yelmox/yelmox_Antarctica_paleo_lgp.nml -o "${output_path}" \
      -p ctrl.restart="${restart}"
