#!/usr/bin/env bash
# Build + run the snapclim validation-reference dumper (tests/test_snapclim_ref.f90).
# Produces logs/snapclim_ref.nc: the tight-tolerance reference for the snapclim2 port.
# Run from the yelmox root (relative paths yelmox/, ice_data/, input/, logs/ must resolve).
set -euo pipefail

INC=fesm-utils/include-serial
BUILD=$(mktemp -d)
NCF=/opt/homebrew/Cellar/netcdf-fortran/4.6.2/lib
NCC=/opt/homebrew/Cellar/netcdf/4.10.0/lib
mkdir -p logs

gfortran -ffree-line-length-none -O2 -I"$INC" -c libs/snapclim.f90 -o "$BUILD/snapclim.o" -J"$BUILD"

gfortran -ffree-line-length-none -O2 -I"$BUILD" -I"$INC" \
    tests/test_snapclim_ref.f90 "$BUILD/snapclim.o" \
    -L"$INC" -lfesmutils \
    -L"$NCF" -lnetcdff -L"$NCC" -lnetcdf \
    -Wl,-rpath,"$NCF" -Wl,-rpath,"$NCC" \
    -o "$BUILD/test_snapclim_ref.x"

"$BUILD/test_snapclim_ref.x"
