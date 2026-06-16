# Compiler configuration: GNU Fortran (gfortran).
#
# A machine fragment (config/machines/<machine>.mk) is loaded *after* this
# file and may override any of these variables (the later assignment wins) --
# e.g. DFLAGS_NODEBUG for CPU-specific optimization.

FC = gfortran

FFLAGS = -ffree-line-length-none -I$(objdir) -J$(objdir)
FFLAGS_OPENMP = -fopenmp

DFLAGS_NODEBUG = -O2
# Trap real FP errors (NaN, div-by-zero, overflow) but NOT underflow: legitimate
# denormal intermediates (e.g. strain-rate sqrt in deformation.f90) otherwise raise
# spurious SIGFPE during normal, stable runs.
DFLAGS_DEBUG   = -w -g -ggdb -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all
DFLAGS_PROFILE = -O2 -pg
