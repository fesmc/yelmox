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
# built yet. They are built one-time by `make libs` inside fesm-utils/, which
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
	  echo "    Build them in fesm-utils with make (openmp=0 serial, openmp=1 omp):"; \
	  echo "      configme config fesm-utils -m <machine> -c <compiler>"; \
	  echo "      make -C fesm-utils libs openmp=0   # then openmp=1 for the omp variant"; \
	  echo "    e.g. -m dkrz_levante -c ifx  (see configme list for machines/compilers)"; \
	  echo ""; \
	  exit 1; \
	fi

# Compile static library in its own directory
fesmutils-static: fesmutils-deps-check
	$(MAKE) -C $(FESMUTILSROOT) fesmutils-static openmp=$(openmp)

# Compile static library in its own directory.
# yelmo's own Makefile builds FastHydrology (yelmo's dependency) as part of
# yelmo-static; fesmutils-static must run first so the fesm-utils/FFTW includes
# FastHydrology needs are in place before yelmo descends into it.
yelmo-static: fesmutils-static
	$(MAKE) -C $(YELMOROOT) yelmo-static openmp=$(openmp)
	cp $(YELMOROOT)/git_yelmo.txt ./

# Compile static library in its own directory
isostasy-static: fesmutils-static
	$(MAKE) -C $(ISOSTASYROOT) isostasy-static openmp=$(openmp)

# Compile the static library librembo,
# using Makefile located in $(REMBOROOT) directory.
# rembo1 builds against fesm-utils (which provides the coord/grid modules), so
# fesmutils-static must run first.
rembo-static: fesmutils-static
	$(MAKE) -C $(REMBOROOT) rembo-static

# ---------------------------------------------------------------------------
# Default programs: the new multigrid drivers, one per flavor folder.
# Legacy single-grid originals live under <flavor>/legacy/ and build via the
# `<name>-legacy` targets below. Flavors not yet ported print a stub.
# ---------------------------------------------------------------------------

yelmox : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox.x yelmox/yelmox.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox.x is ready."
	@echo " "

yelmox_bipolar : fesmutils-static yelmo-static isostasy-static $(yelmox_libs) $(obm_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_bipolar.x yelmox_bipolar/yelmox_bipolar.f90 \
		$(yelmox_libs) $(obm_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_bipolar.x is ready."
	@echo " "

yelmox_esm : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -c -o $(objdir)/yelmox_esm_output.o yelmox_esm/yelmox_esm_output.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_esm.x yelmox_esm/yelmox_esm.f90 \
		$(objdir)/yelmox_esm_output.o $(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_esm.x is ready."
	@echo " "

# REMBO-coupled driver: reuses the multigrid ice_domain, with REMBO as the
# atmosphere/SMB (step_rembo) and its own combined output module.
yelmox_rembo : fesmutils-static yelmo-static rembo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) $(INC_REMBO) -c -o $(objdir)/yelmox_rembo_output.o yelmox_rembo/yelmox_rembo_output.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) $(INC_REMBO) -o $(bindir)/yelmox_rembo.x yelmox_rembo/yelmox_rembo.f90 \
		$(objdir)/yelmox_rembo_output.o $(yelmox_libs) $(LIB_REMBO) $(LFLAGS)
	@echo " "
	@echo "    yelmox_rembo.x is ready."
	@echo " "

# --- flavors not yet ported to the multigrid driver: stub default targets ---
.PHONY : yelmox_ismip6 yelmox_nahosmip yelmox_rtip
yelmox_ismip6 yelmox_nahosmip yelmox_rtip :
	@echo ""
	@echo "    *** $@ is not ported to the new (multigrid) driver yet."
	@echo "    Build the legacy program instead:  make $@-legacy"
	@echo ""

# ---------------------------------------------------------------------------
# Legacy programs (single-grid originals under <flavor>/legacy/).
# ---------------------------------------------------------------------------

yelmox-legacy : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox-legacy.x yelmox/legacy/yelmox.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox-legacy.x is ready."
	@echo " "

yelmox_bipolar-legacy : fesmutils-static yelmo-static isostasy-static $(yelmox_libs) $(obm_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_bipolar-legacy.x yelmox_bipolar/legacy/yelmox_bipolar.f90 \
		$(yelmox_libs) $(obm_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_bipolar-legacy.x is ready."
	@echo " "

yelmox_esm-legacy : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_esm-legacy.x yelmox_esm/legacy/yelmox_esm.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_esm-legacy.x is ready."
	@echo " "

yelmox_ismip6-legacy : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_ismip6-legacy.x yelmox_ismip6/legacy/yelmox_ismip6.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_ismip6-legacy.x is ready."
	@echo " "

yelmox_nahosmip-legacy : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_nahosmip-legacy.x yelmox_nahosmip/legacy/yelmox_nahosmip.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_nahosmip-legacy.x is ready."
	@echo " "

yelmox_rtip-legacy : fesmutils-static yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_rtip-legacy.x yelmox_rtip/legacy/yelmox_rtip.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_rtip-legacy.x is ready."
	@echo " "

yelmox_rembo-legacy : fesmutils-static yelmo-static rembo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) $(INC_REMBO) -o $(bindir)/yelmox_rembo-legacy.x yelmox_rembo/legacy/yelmox_rembo.f90 \
		$(yelmox_libs) $(LIB_REMBO) $(LFLAGS)
	@echo " "
	@echo "    yelmox_rembo-legacy.x is ready."
	@echo " "


test_timestepping : fesmutils-static yelmo-static $(objdir)/nml.o $(objdir)/ncio.o $(objdir)/timestepping.o
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_YELMO) -o ./libyelmox/bin/test_timestepping.x tests/test_timestepping.f90 \
		$(objdir)/nml.o $(objdir)/ncio.o $(objdir)/timestepping.o $(LFLAGS)
	@echo " "
	@echo "    ./libyelmox/bin/test_timestepping.x is ready."
	@echo " "


test_htopo : fesmutils-static yelmo-static isostasy-static $(objdir)/htopo.o
	$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_htopo.x tests/test_htopo.f90 \
		$(objdir)/htopo.o $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_htopo.x is ready."
	@echo " "

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " Default (multigrid) drivers, one per flavor folder:"
	@echo "   make yelmox | yelmox_bipolar | yelmox_esm | yelmox_rembo"
	@echo " Flavors not yet ported (print a stub, build legacy instead):"
	@echo "   make yelmox_ismip6 | yelmox_nahosmip | yelmox_rtip"
	@echo " Legacy single-grid originals (in <flavor>/legacy/):"
	@echo "   make <name>-legacy   e.g. make yelmox-legacy, make yelmox_rembo-legacy"
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
	$(MAKE) -C $(REMBOROOT) clean
