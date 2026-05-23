# Machine configuration: awi_albedo (intended compiler: ifx).
NETCDFC_ROOT = /albedo/soft/sw/spack-sw/netcdf-c/4.8.1-5ewdrxn
NETCDFFI_ROOT = /albedo/soft/sw/spack-sw/netcdf-fortran/4.5.4-rujc52i
INC_NC = -I${NETCDFFI_ROOT}/include
LIB_NC = -L${NETCDFFI_ROOT}/lib -Wl\,-rpath=${NETCDFFI_ROOT}/lib -lnetcdff -L${NETCDFC_ROOT}/lib -Wl\,-rpath=${NETCDFC_ROOT}/lib -lnetcdf
# the -rpath options ensure that the same shared libraries are found at runtime

# CPU-specific optimization; overrides the ifx default DFLAGS_NODEBUG.
DFLAGS_NODEBUG = -Ofast -march=core-avx2 -mtune=core-avx2 -traceback
