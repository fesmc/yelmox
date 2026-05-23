# Machine configuration: anta (intended compiler: gfortran).
# System netCDF under /usr.
INC_NC = -I/usr/include
LIB_NC = -L/usr/lib -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lz -lcurl -lm
