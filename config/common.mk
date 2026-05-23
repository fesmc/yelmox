# Shared build configuration.
#
# Loaded *after* the compiler and machine fragments (config.py assembles them
# in the order: compiler -> machine -> common). This file references variables
# those fragments define: FFLAGS / FFLAGS_OPENMP (compiler), LIB_NC (machine),
# and LFLAGS_EXTRA (compiler, optionally extended by a machine).

# Dependency paths (serial build by default).
FESMUTILSROOT = fesm-utils/utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

FFTWROOT = fesm-utils/fftw-serial
INC_FFTW = -I${FFTWROOT}/include
LIB_FFTW = -L${FFTWROOT}/lib -lfftw3 -lm

LISROOT = fesm-utils/lis-serial
INC_LIS = -I${LISROOT}/include
LIB_LIS = -L${LISROOT}/lib -llis

ISOSTASYROOT = FastIsostasy
INC_ISOSTASY = -I${ISOSTASYROOT}/libisostasy/include
LIB_ISOSTASY = -L${ISOSTASYROOT}/libisostasy/include -lisostasy

YELMOROOT = yelmo
INC_YELMO = -I${YELMOROOT}/libyelmo/include
LIB_YELMO = -L${YELMOROOT}/libyelmo/include -lyelmo

REMBOROOT = rembo1
INC_REMBO = -I${REMBOROOT}/librembo/include
LIB_REMBO = -L${REMBOROOT}/librembo/include -lrembo

# coordinates is a build dependency of rembo1; only used by `make yelmox_rembo`.
COORDROOT = coordinates

# OpenMP build (make openmp=1): swap the serial dependency builds for their
# OpenMP variants and append the compiler's OpenMP flag (set in the compiler
# fragment as FFLAGS_OPENMP).
ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils

    FFTWROOT = fesm-utils/fftw-omp
    INC_FFTW = -I${FFTWROOT}/include
    LIB_FFTW = -L${FFTWROOT}/lib -lfftw3_omp -lfftw3 -lm

    LISROOT = fesm-utils/lis-omp
    INC_LIS = -I${LISROOT}/include
    LIB_LIS = -L${LISROOT}/lib -llis

    FFLAGS += $(FFLAGS_OPENMP)
endif

LFLAGS = $(LIB_YELMO) $(LIB_ISOSTASY) $(LIB_FESMUTILS) $(LIB_NC) $(LIB_LIS) $(LIB_FFTW) $(LFLAGS_EXTRA)
