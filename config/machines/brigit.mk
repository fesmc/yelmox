# Machine configuration: brigit (intended compiler: ifort).
INC_NC = -I/opt/ohpc/pub/libs/intel/impi/netcdf-fortran/4.4.5/include
LIB_NC = -L/opt/ohpc/pub/libs/intel/impi/netcdf-fortran/4.4.5/lib -lnetcdff -L/opt/ohpc/pub/libs/intel/impi/hdf5/1.10.5/lib -L/opt/ohpc/pub/libs/intel/impi/netcdf/4.6.3/lib -lnetcdf -lnetcdf -lhdf5_hl -lhdf5 -lz -lm
# Uses the ifort default DFLAGS_NODEBUG (-O2 -fp-model precise); no CPU override.
