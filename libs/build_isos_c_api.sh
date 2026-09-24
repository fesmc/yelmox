#!/usr/bin/env bash
# Builds libyelmox_isos_c_api.so from libs/yelmox_isos_c_api.f90.
#
# Mirrors build_marshelf_c_api.sh (see its header comment for the general
# rationale: top-level Makefile is configme-generated/gitignored, so this
# is the trackable build path). Run after `configme install yelmox` /
# `make yelmox_esm` has already built yelmox's usual dependencies.
#
# Usage: from the yelmox root (or a worktree of it), with the same
# compiler module loaded that built yelmox_esm.x:
#   module load intel-oneapi-compilers/2024.1.0   # or your machine's equivalent
#   module load fftw/3.3.10                       # for the isostasy FFT solver
#   ./libs/build_isos_c_api.sh
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root

FC=${FC:-ifx}
FFLAGS="-Ofast -march=znver2 -traceback -fPIC -no-wrap-margin"
OBJDIR=libyelmox/include
INC="-module $OBJDIR -L$OBJDIR -Ifesm-utils/include-serial -Ifesm-utils/lis/lis-serial/include -IFastIsostasy/libisostasy/include -Iyelmo/libyelmo/include"

YELMOX_LIBS=$(awk '/^yelmox_libs/{f=1} f{print} /^$/{if(f)exit}' config/Makefile_yelmox.mk \
  | grep -oE '\$\(objdir\)/[A-Za-z0-9_]+\.o' | sed "s#\$(objdir)#$OBJDIR#")

for f in $YELMOX_LIBS; do
  [ -f "$f" ] || { echo "Missing $f -- run 'make yelmox_esm' first." >&2; exit 1; }
done

# Unlike marine_shelf.f90, fastisostasy.f90 genuinely uses FFTW (its
# viscoelastic/convolution solver). fesm-utils' own fftw-serial static lib
# was not built -fPIC (see build_marshelf_c_api.sh's note) so it can't link
# into a .so -- use the system fftw/3.3.10 module's shared libfftw3.so
# instead, which has no such restriction.
FFTW_LIBDIR=$(module show fftw/3.3.10 2>&1 | grep LD_LIBRARY_PATH | awk '{print $3}')
[ -n "${FFTW_LIBDIR:-}" ] && [ -d "$FFTW_LIBDIR" ] || { echo "Could not resolve fftw/3.3.10 lib dir -- 'module load fftw/3.3.10' first?" >&2; exit 1; }

$FC -shared -fPIC $FFLAGS $INC -c -o "$OBJDIR/yelmox_isos_c_api.o" libs/yelmox_isos_c_api.f90

$FC -shared -fPIC $FFLAGS $INC \
  -o "$OBJDIR/libyelmox_isos_c_api.so" \
  "$OBJDIR/yelmox_isos_c_api.o" $YELMOX_LIBS \
  -Lyelmo/libyelmo/include -lyelmo \
  -LFastIsostasy/libisostasy/include -lisostasy \
  -Lfesm-utils/include-serial -lfesmutils \
  -Lyelmo/FastHydrology/include -lfasthydro \
  -Lyelmo/elsa/libelsa/include -lelsa \
  -Lyelmo/tracer/libtracer/include -ltracer \
  -L"$FFTW_LIBDIR" -lfftw3 -Wl,-rpath,"$FFTW_LIBDIR" \
  $(nf-config --flibs 2>/dev/null || echo "-L/albedo/soft/sw/spack-sw/netcdf-fortran/4.5.4-lzqfsg3/lib -lnetcdff -L/albedo/soft/sw/spack-sw/netcdf-c/4.8.1-5ewdrxn/lib -lnetcdf") \
  -Wl,-zmuldefs

echo
echo "    $OBJDIR/libyelmox_isos_c_api.so is ready."
echo
