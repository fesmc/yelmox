#!/usr/bin/env python3
"""Check yelmox parameter files for keys that break or mislead a run.

Yelmo reads each parameter with `nml_read(..., defaults_file=...)`: a
parameter missing from the par file takes its value from
`input/yelmo_defaults.nml`. Before reading a group, `nml_validate` checks every
key the par file sets in that group against the same group of the defaults
file and stops the run on any key it does not know (typos, removed or renamed
parameters). This script reproduces that check offline, with two more:

  1. unknown: keys in a yelmo group that are not in the defaults file
     (`nml_validate` would stop the run).
  2. duplicate: a group or key defined twice in one file. `nml_read` takes the
     first occurrence, so later copies are silently ignored.
  3. dead (only with --dead): keys in non-yelmo groups (ctrl, snapclim,
     smbpal, isos, ...) whose name appears nowhere in the Fortran sources as a
     string literal, so no `nml_read` can read them. The match is by key name
     only (not by group), so it misses keys that are read in another group.
     Needs the dependency checkouts (yelmo, fesm-utils, FastIsostasy, rembo1)
     next to the scripts directory.

The yelmo groups of a file are found from its control blocks: `&yelmo<sfx>`
for every `&domain<sfx>` group (kryos `domain_init` reads `yelmo` plus the
domain's group suffix), plus `&yelmo` itself. The `nml_*` keys of a control
block name its component groups (`nml_ytopo = "ytopo"`, ...), falling back to
the defaults. Files without a control block are skipped.

Usage:
    scripts/check_par_nml.py                  # all yelmox/*.nml, yelmox_bipolar/*.nml
    scripts/check_par_nml.py yelmox/yelmox_Greenland.nml [...]
    scripts/check_par_nml.py --dead           # also report dead keys

Exits non-zero if any error is found, so it can be used as a test / CI check.
"""
import argparse
import glob
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULTS = os.path.join(REPO, "input", "yelmo_defaults.nml")
PAR_GLOBS = ["yelmox/*.nml", "yelmox_bipolar/*.nml"]
# Source trees searched by --dead (symlinked checkouts are followed, so yelmo
# brings its own elsa, tracer, FastHydrology, ...).
SRC_DIRS = ["libs", "yelmox", "yelmox_bipolar", "fesm-utils", "FastIsostasy",
            "rembo1", "yelmo"]
SKIP_DIRS = {".git", ".claude", "output", "logs", "tmp"}

# nml_* key of the control block -> group of the defaults file it maps to.
NML_KEY = {
    "nml_ytopo": "ytopo", "nml_ycalv": "ycalv", "nml_ydyn": "ydyn",
    "nml_ytill": "ytill", "nml_ymat": "ymat", "nml_ytrc": "ytrc",
    "nml_ytherm": "ytherm", "nml_yhyd": "yhyd", "nml_masks": "yelmo_masks",
    "nml_init_topo": "yelmo_init_topo", "nml_data": "yelmo_data",
}


def parse_par(path):
    """Parse a par file like fesm-utils nml.f90 (case-sensitive names).

    Returns (groups, errors): groups is a list of (name, {key: value}) in file
    order, keeping duplicate groups; errors lists duplicate keys and lines
    nml_read cannot parse.
    """
    groups, errors = [], []
    cur = None
    for i, line in enumerate(open(path), 1):
        s = line.strip()
        if not s or s.startswith("!"):
            continue
        if s.startswith("&"):
            cur = {}
            groups.append((s[1:].strip(), cur))
            continue
        if s.startswith("/"):
            cur = None
            continue
        if "=" not in s:
            errors.append("line %d: no '=' in parameter line" % i)
            continue
        if cur is None:
            continue
        key, val = s.split("=", 1)
        key = key.strip()
        if key in cur:
            errors.append("&%s: duplicate key %s (line %d ignored)"
                          % (groups[-1][0], key, i))
            continue
        cur[key] = val.split("!", 1)[0].strip().rstrip(",").strip("'\"")
    return groups, errors


def first_copies(groups):
    """Group name -> keys of its first copy, plus errors for later copies."""
    out, errors = {}, []
    for name, params in groups:
        if name in out:
            errors.append("&%s: duplicate group (later copy ignored)" % name)
        else:
            out[name] = params
    return out, errors


def control_blocks(groups):
    """Names of the &yelmo<sfx> groups kryos passes to yelmo_init."""
    sfxs = {""} | {g[len("domain"):] for g in groups if g.startswith("domain")}
    return [("yelmo" + s) for s in sorted(sfxs) if ("yelmo" + s) in groups]


def yelmo_groups(groups, defaults):
    """User group name -> defaults group name, for all yelmo groups of a file."""
    out = {}
    for ctrl in control_blocks(groups):
        out[ctrl] = "yelmo"
        for key, def_group in NML_KEY.items():
            name = groups[ctrl].get(key, defaults["yelmo"].get(key, def_group))
            out.setdefault(name, def_group)
    return out


def source_literals():
    """All identifier-like string literals in the Fortran sources."""
    lits = set()
    for d in SRC_DIRS:
        for root, dirs, files in os.walk(os.path.join(REPO, d), followlinks=True):
            dirs[:] = [x for x in dirs if x not in SKIP_DIRS]
            for f in files:
                if not f.lower().endswith((".f90", ".f")):
                    continue
                for line in open(os.path.join(root, f), errors="replace"):
                    if line.lstrip().startswith("!"):
                        continue
                    lits.update(m[0] or m[1] for m in re.findall(
                        r'"([A-Za-z_]\w*)"|\'([A-Za-z_]\w*)\'', line))
    return lits


def check_file(path, defaults, lits):
    """Return (errors, dead) for one par file, or None if it has no yelmo group."""
    raw, errors = parse_par(path)
    groups, dup = first_copies(raw)
    if not control_blocks(groups):
        return None
    errors += dup

    ygroups = yelmo_groups(groups, defaults)
    for name, def_group in ygroups.items():
        unknown = [k for k in groups.get(name, {}) if k not in defaults[def_group]]
        if unknown:
            errors.append("&%s: unknown key(s) (not in &%s of defaults): %s"
                          % (name, def_group, ", ".join(unknown)))

    dead = []
    if lits is not None:
        for name, params in groups.items():
            if name in ygroups:
                continue
            keys = [k for k in params if k not in lits]
            if keys:
                dead.append("&%s: %s" % (name, ", ".join(keys)))
    return errors, dead


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*",
                    help="par files to check (default: %s)" % ", ".join(PAR_GLOBS))
    ap.add_argument("--defaults", default=DEFAULTS,
                    help="yelmo defaults file (default: input/yelmo_defaults.nml)")
    ap.add_argument("--dead", action="store_true",
                    help="also report non-yelmo keys no source reads (errors)")
    args = ap.parse_args(argv[1:])

    if not os.path.isfile(args.defaults):
        print("ERROR: defaults file not found: %s" % args.defaults)
        return 2
    defaults, def_errors = parse_par(args.defaults)
    defaults, dup = first_copies(defaults)
    if def_errors or dup:
        print("ERROR: defaults file: " + "; ".join(def_errors + dup))
        return 2

    lits = None
    if args.dead:
        missing = [d for d in SRC_DIRS if not os.path.isdir(os.path.join(REPO, d))]
        if missing:
            print("ERROR: --dead needs the source trees: %s" % ", ".join(missing))
            return 2
        lits = source_literals()
    files = args.files or sorted(f for g in PAR_GLOBS
                                 for f in glob.glob(os.path.join(REPO, g)))

    print("Defaults: %s" % os.path.relpath(args.defaults, REPO))
    print()

    n_fail = 0
    for path in files:
        rel = os.path.relpath(path, REPO)
        result = check_file(path, defaults, lits)
        if result is None:
            print("[skip] %-40s no &yelmo control block" % rel)
            continue
        errors, dead = result
        if not errors and not dead:
            print("[ ok ] %s" % rel)
            continue
        n_fail += 1
        print("[FAIL] %s" % rel)
        for msg in errors:
            print("         " + msg)
        for msg in dead:
            print("         dead " + msg)

    print()
    if n_fail:
        print("%d file(s) with errors." % n_fail)
        return 1
    print("All par files OK.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
