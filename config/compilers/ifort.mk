# Compiler configuration: Intel Fortran classic (ifort).
#
# A machine fragment (config/machines/<machine>.mk) is loaded *after* this
# file and may override any of these variables (the later assignment wins).
# DFLAGS_NODEBUG here is a safe, portable default; machines with a known CPU
# may override it with -march tuning.

FC = ifort

FFLAGS = -no-wrap-margin -module $(objdir) -L$(objdir)
FFLAGS_OPENMP = -qopenmp

DFLAGS_NODEBUG = -O2 -fp-model precise
DFLAGS_DEBUG   = -C -O0 -g -traceback -ftrapuv -fpe0 -check all,nouninit -fp-model precise -debug extended -gen-interfaces -warn interfaces -check arg_temp_created
DFLAGS_PROFILE = -O2 -fp-model precise -pg
