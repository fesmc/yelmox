# Machine configuration: dkrz_levante (intended compiler: ifx).
NETCDFC_ROOT = /sw/spack-levante/netcdf-c-4.8.1-2k3cmu
NETCDFFI_ROOT = /sw/spack-levante/netcdf-fortran-4.5.3-k6xq5g
INC_NC = -I${NETCDFFI_ROOT}/include
LIB_NC = -L${NETCDFFI_ROOT}/lib -Wl\,-rpath=${NETCDFFI_ROOT}/lib -lnetcdff -L${NETCDFC_ROOT}/lib -Wl\,-rpath=${NETCDFC_ROOT}/lib -lnetcdf
# the -rpath options ensure that the same shared libraries are found at runtime

# CPU-specific optimization; overrides the ifx default DFLAGS_NODEBUG.
DFLAGS_NODEBUG = -Ofast -march=core-avx2 -mtune=core-avx2 -traceback
