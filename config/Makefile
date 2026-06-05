.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = libyelmox/include
bindir = libyelmox/bin
libdir = libs

# Command-line options at make call
debug  ?= 0
openmp ?= 0 

## COMPILER CONFIGURATION ##
# Loaded from a config/<host> fragment by `python config.py --file <host>`,
# or assembled by `python config.py <machine> <compiler>` (modern path).
# The fragment / pair defines the compiler (FC, FFLAGS, FFLAGS_OPENMP,
# DFLAGS_*), netCDF wiring (INC_NC, LIB_NC), and any host-specific
# overrides (e.g. LFLAGS_EXTRA).

<COMPILER_CONFIGURATION>

## SHARED DEPENDENCY WIRING ##
# Pulled in after the compiler/machine fragment so it can reference
# FFLAGS, FFLAGS_OPENMP, LIB_NC, and LFLAGS_EXTRA set there. Defines the
# yelmo, FastIsostasy, FESMUTILS, LIS, FastHydrology, FFTW, and rembo1
# paths, the openmp= build switch, and the final LFLAGS.

include config/common.mk

# Determine whether to use normal flags or debugging flags
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
endif

###############################################
##
## List of yelmox rules and source files
##
###############################################

include config/Makefile_yelmox.mk

###############################################
##
## Compilation of complete programs
##
###############################################

# Sanity check: fail early if fesm-utils' external deps (LIS + FFTW) aren't
# built yet. They are built one-time by ./build.py inside fesm-utils/, which
# is kept as a manual step because compiling LIS + FFTW is slow (~10-30 min)
# and machine/compiler-specific.
fesmutils-deps-check:
	@_lis_ok=1; _fftw_ok=1; \
	[ -f $(LISROOT)/lib/liblis.a ]   || _lis_ok=0; \
	[ -f $(FFTWROOT)/lib/libfftw3.a ] || _fftw_ok=0; \
	if [ $$_lis_ok -eq 0 ] || [ $$_fftw_ok -eq 0 ]; then \
	  echo ""; \
	  echo "*** ERROR: fesm-utils external deps not built:"; \
	  [ $$_lis_ok  -eq 0 ] && echo "      LIS  missing at $(LISROOT)/lib/liblis.a"; \
	  [ $$_fftw_ok -eq 0 ] && echo "      FFTW missing at $(FFTWROOT)/lib/libfftw3.a"; \
	  echo ""; \
	  echo "    Build them in fesm-utils with build.py (--variant both builds serial + omp):"; \
	  echo "      ./build.py --variant both --machine <machine> --compiler <compiler>"; \
	  echo "    Supported machines + compilers (./build.py --list-machines for the live list):"; \
	  echo "      awi_albedo, dkrz_levante, generic, pik_hpc2024 : ifx, ifort, gfortran"; \
	  echo "      macbook                                        : gfortran"; \
	  echo "    e.g. ./build.py --variant both --machine dkrz_levante --compiler ifx"; \
	  echo ""; \
	  exit 1; \
	fi

# Compile static library in its own directory
fesmutils-static: fesmutils-deps-check
	$(MAKE) -C $(FESMUTILSROOT) fesmutils-static openmp=$(openmp)

# Compile static library in its own directory
yelmo-static: fesmutils-static fasthydro-static
	$(MAKE) -C $(YELMOROOT) yelmo-static openmp=$(openmp)
	cp $(YELMOROOT)/git_yelmo.txt ./

# FastHydrology's static archive provides fast_hydrology.mod, which yelmo's
# sources `use`, so this must build before yelmo-static.
fasthydro-static: fesmutils-static
	$(MAKE) -C $(FASTHYDROROOT) fasthydro-static openmp=$(openmp)

# Compile static library in its own directory
isostasy-static: fesmutils-static
	$(MAKE) -C $(ISOSTASYROOT) isostasy-static openmp=$(openmp)

# Compile static library in its own directory
coord-static:
	$(MAKE) -C $(COORDROOT) coord-static

# Compile the static library librembo,
# using Makefile located in $(REMBOROOT) directory
rembo-static: coord-static
	$(MAKE) -C $(REMBOROOT) rembo-static

yelmox : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox.x yelmox.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox.x is ready."
	@echo " "

yelmoxmg : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmoxmg.x yelmox.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmoxmg.x is ready."
	@echo " "

yelmox_hyst : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_hyst.x yelmox_hyst.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_hyst.x is ready."
	@echo " "

yelmox_iso : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_iso.x yelmox_iso.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_iso.x is ready."
	@echo " "

yelmox_ismip6 : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_ismip6.x yelmox_ismip6.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_ismip6.x is ready."
	@echo " "

yelmox_esm : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_esm.x yelmox_esm.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_esm.x is ready."
	@echo " "

yelmox_glaciers : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_glaciers.x yelmox_glaciers.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_glaciers.x is ready."
	@echo " "

yelmox_rtip : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_rtip.x yelmox_rtip.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_rtip.x is ready."
	@echo " "

yelmox_nahosmip : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_nahosmip.x yelmox_nahosmip.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_nahosmip.x is ready."
	@echo " "

yelmox_rembo: fesmutils-static yelmo-static rembo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) $(INC_REMBO) -o libyelmox/bin/yelmox_rembo.x \
			yelmox_rembo.f90 $(LIB_REMBO) $(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    libyelmox/bin/yelmox_rembo.x is ready."
	@echo " "

yelmox_bipolar : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_YELMO) $(INC_ISOSTASY) -o $(bindir)/yelmox_bipolar.x yelmox_bipolar.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_bipolar.x is ready."
	@echo " "

check_sim : $(objdir)/ncio.o
	$(FC) $(DFLAGS) $(FFLAGS) -o ./check_sim.x check_sim.f90 \
		$(objdir)/ncio.o $(LFLAGS)
	@echo " "
	@echo "    check_sim.x is ready."
	@echo " "

test_timestepping : fesmutils-static yelmo-static $(objdir)/nml.o $(objdir)/ncio.o $(objdir)/timestepping.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_YELMO) -o ./libyelmox/bin/test_timestepping.x tests/test_timestepping.f90 \
		$(objdir)/nml.o $(objdir)/ncio.o $(objdir)/timestepping.o $(LFLAGS)
	@echo " "
	@echo "    ./libyelmox/bin/test_timestepping.x is ready."
	@echo " "


.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make yelmox     : compiles yelmox.x, for running yelmo on a given domain defined in param file."
	@echo " make clean      : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM

cleanall: clean
	$(MAKE) -C $(FESMUTILSROOT) clean
	$(MAKE) -C $(ISOSTASYROOT) clean
	$(MAKE) -C $(YELMOROOT) clean
	$(MAKE) -C $(COORDROOT) clean
	$(MAKE) -C $(REMBOROOT) clean
