# Machine configuration: macbook (intended compiler: gfortran).
# NetCDF roots come from the environment; set NC_CROOT and NC_FROOT in .zshrc.

INC_NC = -I${NC_FROOT}/include
LIB_NC = -L${NC_FROOT}/lib -lnetcdff -L${NC_CROOT}/lib -lnetcdf

# Disable the default -Wl,-zmuldefs: gfortran forwards it to Apple's ld (ld64),
# which rejects it ("ld: unknown options: -zmuldefs"). It is a GNU-ld/ELF flag.
LFLAGS_EXTRA =

