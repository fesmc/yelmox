# Machine configuration: pik_hpc2024 (intended compiler: ifx).
# NETCDFC_ROOT / NETCDFFI_ROOT are provided by the environment (module load).
INC_NC = -I${NETCDFFI_ROOT}/include
LIB_NC = -L${NETCDFFI_ROOT}/lib -Wl\,-rpath=${NETCDFFI_ROOT}/lib -lnetcdff -L${NETCDFC_ROOT}/lib -Wl\,-rpath=${NETCDFC_ROOT}/lib -lnetcdf
# the -rpath options ensure that the same shared libraries are found at runtime
# Uses the ifx default DFLAGS_NODEBUG (-O2 -fp-model precise); no CPU override.
