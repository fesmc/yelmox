# Machine configuration: macbook (intended compiler: gfortran).
# NetCDF roots come from the environment; set NC_CROOT and NC_FROOT in .zshrc.
INC_NC = -I${NC_FROOT}/include
LIB_NC = -L${NC_FROOT}/lib -lnetcdff -L${NC_CROOT}/lib -lnetcdf

# macOS ld does not support -Wl,-zmuldefs; clear the compiler default.
LFLAGS_EXTRA =
