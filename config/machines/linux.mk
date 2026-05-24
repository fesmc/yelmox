# Machine configuration: linux (intended compiler: gfortran).
# NetCDF roots come from the environment; set NC_CROOT and NC_FROOT in .bashrc.

INC_NC = -I${NC_FROOT}/include
LIB_NC = -L${NC_FROOT}/lib -lnetcdff -L${NC_CROOT}/lib -lnetcdf

