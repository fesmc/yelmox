# Compiler configuration: GNU Fortran (gfortran).
#
# A machine fragment (config/machines/<machine>.mk) is loaded *after* this
# file and may override any of these variables (the later assignment wins) --
# e.g. DFLAGS_NODEBUG for CPU-specific optimization.

FC = gfortran

FFLAGS = -ffree-line-length-none -I$(objdir) -J$(objdir)
FFLAGS_OPENMP = -fopenmp

# -Wl,-zmuldefs works around duplicate symbols when linking the static deps.
# Needed on the Linux target machines; macOS ld lacks it, so macbook overrides
# this to empty in its machine fragment.
LFLAGS_EXTRA = -Wl,-zmuldefs

DFLAGS_NODEBUG = -O2
DFLAGS_DEBUG   = -w -g -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
DFLAGS_PROFILE = -O2 -pg
