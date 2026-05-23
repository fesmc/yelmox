# Machine configuration: foehn (intended compiler: gfortran).
# System netCDF under /usr.
INC_NC = -I/usr/include
LIB_NC = -L/usr/lib/x86_64-linux-gnu -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lm
