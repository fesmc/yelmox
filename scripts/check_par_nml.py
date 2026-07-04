#!/usr/bin/env python3
"""Check that yelmox parameter files define every yelmo-core namelist parameter.

Yelmo reads its parameters with `nml_read`, which aborts the run
(`ERROR_NO_PARAM = .TRUE.`) if a requested parameter is absent from the
namelist file. A par file missing even one parameter therefore crashes at
startup. This script compares each par file against a reference par file that
holds the complete set of yelmo-core parameters and reports any that are
missing.

The reference defaults to `yelmo/par/yelmo_initmip.nml` (reached via the
`yelmo` symlink at the repo root), which is maintained alongside the yelmo
source. Override it with `--ref`. Only the yelmo-core groups present in the
reference are checked; flavor/library groups (snapclim, smbpal, marine_shelf,
ismip6, esm, ...) are out of scope.

Usage:
    scripts/check_par_nml.py                      # check all par/*.nml
    scripts/check_par_nml.py yelmox/legacy/yelmo_Greenland.nml [...]
    scripts/check_par_nml.py --ref par/other.nml  # use a different reference

Exits non-zero if any required parameter is missing, so it can be used as a
test / CI check.
"""
import argparse
import glob
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_REF = os.path.join(REPO, "yelmo", "par", "yelmo_initmip.nml")

# Standard group name -> the nml_* key in the control block that remaps it.
# The control block is the &yelmo / &yelmo_* section that defines nml_ytopo etc.
NML_KEY = {
    "ytopo": "nml_ytopo", "ycalv": "nml_ycalv", "ydyn": "nml_ydyn",
    "ytill": "nml_ytill", "yneff": "nml_yneff", "ymat": "nml_ymat",
    "ytherm": "nml_ytherm", "yelmo_data": "nml_data",
}


def parse_par(path):
    """path -> dict group_name -> {param: raw_value_str}."""
    groups = {}
    cur = None
    for line in open(path):
        s = line.strip()
        if not s or s.startswith("!"):
            continue
        if s.startswith("&"):
            cur = s[1:].split()[0]
            groups.setdefault(cur, {})
            continue
        if s.startswith("/"):
            cur = None
            continue
        if cur is None:
            continue
        m = re.match(r'([A-Za-z_]\w*)\s*=\s*(.*)', s)
        if m:
            groups[cur][m.group(1)] = m.group(2)
    return groups


def control_blocks(groups):
    """Control blocks are groups that define nml_ytopo (i.e. &yelmo / &yelmo_*)."""
    return {g: params for g, params in groups.items() if "nml_ytopo" in params}


def unquote(v):
    return v.strip().strip("'\"")


def required_params(ref_path):
    """Build the required parameter set from the reference par file.

    Returns standard group name -> set of parameter names. The "yelmo" key
    holds the control-block (meta) parameters; the remaining keys hold the
    parameters of each core group as named in the reference.
    """
    groups = parse_par(ref_path)
    ctrls = control_blocks(groups)
    if not ctrls:
        raise RuntimeError("reference has no &yelmo control block: " + ref_path)
    # Use the first control block as the reference definition.
    ctrl_params = next(iter(ctrls.values()))

    req = {"yelmo": set(ctrl_params)}
    for std, key in NML_KEY.items():
        actual = unquote(ctrl_params[key]) if key in ctrl_params else std
        req[std] = set(groups.get(actual, {}))
    return req


def check_file(path, req):
    """Return list of (label, sorted_missing_params), or None if not a yelmo par file."""
    groups = parse_par(path)
    ctrls = control_blocks(groups)
    if not ctrls:
        return None

    missing = []
    checked_groups = set()  # avoid double-reporting shared groups

    for ctrl_name, ctrl_params in ctrls.items():
        # The control block itself must contain the "yelmo" meta parameters.
        if ctrl_name not in checked_groups:
            checked_groups.add(ctrl_name)
            miss = sorted(req["yelmo"] - set(ctrl_params))
            if miss:
                missing.append(("&" + ctrl_name, miss))

        # Resolve and check each remapped core group.
        for std, key in NML_KEY.items():
            actual = unquote(ctrl_params[key]) if key in ctrl_params else std
            if actual in checked_groups:
                continue
            checked_groups.add(actual)
            if actual not in groups:
                missing.append(("&%s (group absent)" % actual, sorted(req[std])))
            else:
                miss = sorted(req[std] - set(groups[actual]))
                if miss:
                    missing.append(("&" + actual, miss))
    return missing


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*",
                    help="par files to check (default: all par/*.nml)")
    ap.add_argument("--ref", default=DEFAULT_REF,
                    help="reference par file defining the required parameters "
                         "(default: yelmo/par/yelmo_initmip.nml)")
    args = ap.parse_args(argv[1:])

    if not os.path.isfile(args.ref):
        print("ERROR: reference par file not found: %s" % args.ref)
        return 2

    req = required_params(args.ref)
    files = args.files or sorted(glob.glob(os.path.join(REPO, "par", "*.nml")))

    print("Reference: %s" % os.path.relpath(args.ref, REPO))
    print()

    n_fail = 0
    for path in files:
        rel = os.path.relpath(path, REPO)
        result = check_file(path, req)
        if result is None:
            print("[skip] %-34s not a yelmo par file" % rel)
            continue
        if not result:
            print("[ ok ] %s" % rel)
            continue
        n_fail += 1
        print("[FAIL] %s" % rel)
        for label, miss in result:
            print("         %s: missing %d -> %s" % (label, len(miss), ", ".join(miss)))

    print()
    if n_fail:
        print("%d file(s) missing required yelmo-core parameters." % n_fail)
        return 1
    print("All par files contain the required yelmo-core parameters.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
