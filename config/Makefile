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
# (should be loaded from config directory)

<COMPILER_CONFIGURATION>

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

# Compile static library in its own directory
fesmutils-static: 
	$(MAKE) -C $(FESMUTILSROOT) fesmutils-static openmp=$(openmp)

# Compile static library in its own directory
yelmo-static: 
	$(MAKE) -C $(YELMOROOT) yelmo-static openmp=$(openmp)
	cp $(YELMOROOT)/git_yelmo.txt ./

# Compile static library in its own directory
isostasy-static: 
	$(MAKE) -C $(ISOSTASYROOT) isostasy-static openmp=$(openmp)

# Compile the static library libyelmo,
# using Makefile located in $(YELMOROOT) directory 
rembo-static: 
	$(MAKE) -C $(REMBOROOT) rembo-static 

yelmox : yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox.x yelmox.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox.x is ready."
	@echo " "

yelmoxmg : yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmoxmg.x yelmox.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmoxmg.x is ready."
	@echo " "

yelmox_hyst : yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_hyst.x yelmox_hyst.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_hyst.x is ready."
	@echo " "

yelmox_iso : yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_iso.x yelmox_iso.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_iso.x is ready."
	@echo " "

yelmox_ismip6 : yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_ismip6.x yelmox_ismip6.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_ismip6.x is ready."
	@echo " "

yelmox_esm : yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_esm.x yelmox_esm.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_esm.x is ready."
	@echo " "

yelmox_glaciers : yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_glaciers.x yelmox_glaciers.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_glaciers.x is ready."
	@echo " "

yelmox_rtip : yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_rtip.x yelmox_rtip.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_rtip.x is ready."
	@echo " "

yelmox_nahosmip : yelmo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) -o $(bindir)/yelmox_nahosmip.x yelmox_nahosmip.f90 \
		$(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    yelmox_nahosmip.x is ready."
	@echo " "

yelmox_rembo: yelmo-static rembo-static isostasy-static $(yelmox_libs)
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_FESMUTILS) $(INC_LIS) $(INC_ISOSTASY) $(INC_YELMO) $(INC_REMBO) -o libyelmox/bin/yelmox_rembo.x \
			yelmox_rembo.f90 $(LIB_REMBO) $(yelmox_libs) $(LFLAGS)
	@echo " "
	@echo "    libyelmox/bin/yelmox_rembo.x is ready."
	@echo " "

yelmox_bipolar : yelmo-static isostasy-static $(yelmox_libs)
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

test_timestepping : yelmo-static $(objdir)/nml.o $(objdir)/ncio.o $(objdir)/timestepping.o
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
	$(MAKE) -C $(ISOSTASYROOT) clean
	$(MAKE) -C $(YELMOROOT) clean
	$(MAKE) -C $(REMBOROOT) clean
