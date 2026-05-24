# Compiler configuration: GNU Fortran (gfortran).
#
# A machine fragment (config/machines/<machine>.mk) is loaded *after* this
# file and may override any of these variables (the later assignment wins) --
# e.g. DFLAGS_NODEBUG for CPU-specific optimization.

FC = gfortran

FFLAGS = -ffree-line-length-none -I$(objdir) -J$(objdir)
FFLAGS_OPENMP = -fopenmp

DFLAGS_NODEBUG = -O2
DFLAGS_DEBUG   = -w -g -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
DFLAGS_PROFILE = -O2 -pg
