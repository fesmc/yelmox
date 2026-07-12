#!/usr/bin/env python3
"""Compare two NetCDF files field-by-field for validating a physics port.

Reports per-field max absolute and max relative difference. Also has a
--summary mode to eyeball midpoint values in a single file.
"""
import argparse
import sys

import numpy as np
from netCDF4 import Dataset

# Coordinate/dimension variables that are skipped unless they differ.
COORD_VARS = {"xc", "yc", "x", "y", "depth", "lev", "month", "time"}
# Fields to print per-time midpoint values for in --summary mode.
SUMMARY_FIELDS = ("tas", "pr", "ta_ann", "pr_ann", "to_ann")


def is_float(var):
    return np.issubdtype(var.dtype, np.floating)


# Fill/missing sentinels are excluded from comparison. All physical fields here
# (temperatures in K, precip, elevation, salinity) are strictly positive and below
# 1e4. The floor of 0 also drops partial-missing ocean cells whose 12-month mean is
# pulled negative by a -9999/-9725.85 month (the two codes' missing conventions
# differ, so such contaminated averages are not meaningful to compare).
MISSING_LO = 0.0
MISSING_HI = 1.0e4


def diff_arrays(a, b, atol, mask_missing=True):
    """Return (maxabs, maxrel, n_nonfinite, n_masked) over jointly-valid entries."""
    a = np.ma.filled(np.ma.asarray(a, dtype=np.float64), np.nan)
    b = np.ma.filled(np.ma.asarray(b, dtype=np.float64), np.nan)
    finite = np.isfinite(a) & np.isfinite(b)
    n_nonfinite = int(np.sum(np.isfinite(a) != np.isfinite(b)))
    n_masked = 0
    if mask_missing:
        def phys(x):
            return (x > MISSING_LO) & (x < MISSING_HI)
        valid = finite & phys(a) & phys(b)
        # Count sentinel cells that would otherwise be compared as finite.
        n_masked = int(np.sum(finite & ~valid))
    else:
        valid = finite
    if not np.any(valid):
        return 0.0, 0.0, n_nonfinite, n_masked
    af, bf = a[valid], b[valid]
    absdiff = np.abs(af - bf)
    maxabs = float(np.max(absdiff))
    maxrel = float(np.max(absdiff / (np.abs(bf) + atol)))
    return maxabs, maxrel, n_nonfinite, n_masked


def compare(file_a, file_b, only_vars, rtol, atol, mask_missing=True):
    da, db = Dataset(file_a), Dataset(file_b)
    if only_vars:
        names = only_vars
    else:
        common = set(da.variables) & set(db.variables)
        names = sorted(
            n for n in common if da.variables[n].ndim >= 1 and is_float(da.variables[n])
        )
    n_cmp = n_pass = n_fail = 0
    for name in (only_vars if only_vars else sorted(set(da.variables) | set(db.variables))):
        in_a, in_b = name in da.variables, name in db.variables
        if only_vars and name not in names:
            continue
        if not in_a:
            print(f"{name}  MISSING in A")
            continue
        if not in_b:
            print(f"{name}  MISSING in B")
            continue
        if not only_vars and name not in names:
            continue
        va, vb = da.variables[name], db.variables[name]
        a, b = va[...], vb[...]
        if a.shape != b.shape:
            print(f"{name}  SHAPE MISMATCH a={a.shape} b={b.shape}")
            continue
        # Skip identical coordinate vars unless they differ.
        maxabs, maxrel, n_nonfinite, n_masked = diff_arrays(a, b, atol, mask_missing)
        if name in COORD_VARS and maxabs == 0.0 and n_nonfinite == 0:
            continue
        n_cmp += 1
        ok = not (maxabs > atol and maxrel > rtol)
        status = "PASS" if ok else "FAIL"
        n_pass += ok
        n_fail += not ok
        extra = f"  nonfinite_mismatch={n_nonfinite}" if n_nonfinite else ""
        extra += f"  masked_missing={n_masked}" if n_masked else ""
        print(
            f"{name}  shape={a.shape}  maxabs={maxabs:.6g}  "
            f"maxrel={maxrel:.6g}  {status}{extra}"
        )
    da.close()
    db.close()
    print(f"\nSummary: compared={n_cmp}  PASS={n_pass}  FAIL={n_fail}")
    return 0 if n_fail == 0 else 1


def summarize(path):
    ds = Dataset(path)
    dims = ds.dimensions
    tname = next((d for d in ("time", "month") if d in dims), None)
    for name in sorted(ds.variables):
        if name in COORD_VARS:
            continue
        var = ds.variables[name]
        if var.ndim == 0:
            continue
        arr = var[...]
        # Identify spatial axes (x is last, y is second-to-last by convention).
        has_time = var.dimensions[0] == tname if tname else False
        # midpoint indices on the two trailing (spatial) axes
        idx = [slice(None)] * var.ndim
        for ax in range(var.ndim):
            if var.dimensions[ax] in (tname,):
                continue
            idx[ax] = arr.shape[ax] // 2
        print(f"{name}  shape={arr.shape}  dims={var.dimensions}")
        if name in SUMMARY_FIELDS and has_time:
            nt = arr.shape[0]
            for t in range(nt):
                sidx = list(idx)
                sidx[0] = t
                val = np.ma.filled(arr[tuple(sidx)], np.nan)
                print(f"    t={t}: {float(val):.6g}")
        else:
            val = np.ma.filled(arr[tuple(idx)], np.nan)
            try:
                print(f"    midpoint: {float(val):.6g}")
            except (TypeError, ValueError):
                print(f"    midpoint: {np.asarray(val).ravel()[:5]}")
    ds.close()
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description="Compare two NetCDF files field-by-field.")
    p.add_argument("files", nargs="*", help="file_a.nc file_b.nc")
    p.add_argument("--vars", help="comma-separated list of variables to compare")
    p.add_argument("--rtol", type=float, default=1e-4)
    p.add_argument("--atol", type=float, default=1e-3)
    p.add_argument("--summary", metavar="FILE.NC", help="single-file eyeball mode")
    p.add_argument("--no-mask-missing", action="store_true",
                   help="compare fill/missing sentinel cells too (default: excluded)")
    args = p.parse_args(argv)

    if args.summary:
        if args.files:
            p.error("--summary is mutually exclusive with two-file compare")
        return summarize(args.summary)

    if len(args.files) != 2:
        p.error("provide exactly two files, or use --summary <file.nc>")
    only_vars = args.vars.split(",") if args.vars else None
    return compare(args.files[0], args.files[1], only_vars, args.rtol, args.atol,
                   mask_missing=not args.no_mask_missing)


if __name__ == "__main__":
    sys.exit(main())
