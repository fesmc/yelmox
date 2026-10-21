#!/usr/bin/env python3
"""Check yelmox parameter files for keys that break or mislead a run.

Yelmo reads each parameter with `nml_read(..., defaults_file=...)`: a
parameter missing from the par file takes its value from
`input/yelmo_defaults.nml`. Before reading a group, `nml_validate` checks every
key the par file sets in that group against the same group of the defaults
file and stops the run on any key it does not know (typos, removed or renamed
parameters). The yelmox forcing groups are read the same way against
`input/yelmox_defaults.nml`, named groups against a template group there. This
script reproduces that check offline, with more:

  1. unknown: keys in a yelmo or yelmox forcing group that are not in its
     group of the defaults file (`nml_validate` would stop the run).
  2. duplicate: a group or key defined twice in one file. `nml_read` takes the
     first occurrence, so later copies are silently ignored.
  3. var_defs: a snapshot (or an active [bias], [taper] or [ice_bnd]) binds a
     varslice group that is in neither the par file nor the [forcing] var_defs
     files (the par file is searched first); such groups in the par file are
     not reported by --dead. The series, term and snapshot groups
     [series_<name>], [term_<name>], [snapshot_<name>] are looked up the same
     way: a listed group found in neither is an error, one found in a var_defs
     file is checked there.
  4. dead (only with --dead): keys in the other groups (ctrl, smbpal, isos,
     ...) whose name appears nowhere in the Fortran sources as a string
     literal, so no `nml_read` can read them. The match is by key name
     only (not by group), so it misses keys that are read in another group.
     Needs the dependency checkouts (yelmo,
     fesm-utils, FastIsostasy, rembo1, chion) next to the scripts directory.

The yelmo groups of a file are found from its control blocks: `&yelmo<sfx>`
for every `&domain<sfx>` group (kryos `domain_init` reads `yelmo` plus the
domain's group suffix), plus `&yelmo` itself. The `nml_*` keys of a control
block name its component groups (`nml_ytopo = "ytopo"`, ...), falling back to
the defaults. Files without a control block are skipped, except library
files (see below).

The yelmox forcing groups, with <sfx> the same domain suffixes, and their
group in yelmox_defaults.nml:
  [forcing<sfx>] -> forcing, [series_<name>] -> series, [term_<name>] -> term,
  [snapshots<sfx>] -> snapshots, [snapshot_<name>] -> snapshot, [obm_forcing] -> obm_forcing,
  [taper<sfx>], [bias<sfx>], [orbit<sfx>] -> taper, bias, orbit, [ice_bnd<sfx>] -> ice_bnd;
  the series of [forcing] series, the terms of [forcing<sfx>] terms, the
  reference and the snap0/snap1 of the snapshot terms (with [comps<sfx>]
  climate = snapshots or rembo); unset lists from the defaults.

A library file has no control block and holds only named groups
(snapshot_*, series_*, term_*, elsa_*), shared by par files through the
[forcing] var_defs list. Each group is checked against its template by prefix:
snapshot_* -> snapshot, series_* -> series, term_* -> term; elsa_*
groups (read by elsa) are skipped. Duplicates are reported as in item 2.

Usage:
    scripts/check_par_nml.py                  # all yelmox/*.nml, yelmox_bipolar/*.nml
    scripts/check_par_nml.py yelmox/yelmox_Greenland.nml [...]
    scripts/check_par_nml.py input/forcing_db.nml   # a library file
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
YX_DEFAULTS = os.path.join(REPO, "input", "yelmox_defaults.nml")
PAR_GLOBS = ["yelmox/*.nml", "yelmox_bipolar/*.nml"]
# Source trees searched by --dead (symlinked checkouts are followed, so yelmo
# brings its own elsa, tracer, FastHydrology, ...).
SRC_DIRS = ["libs", "yelmox", "yelmox_bipolar", "fesm-utils", "FastIsostasy",
            "rembo1", "yelmo", "chion"]
SKIP_DIRS = {".git", ".claude", "output", "logs", "tmp"}

# Group prefixes of a library file -> template group of yelmox_defaults.nml
# (None: not checked here). First match wins.
LIBRARY_PREFIX = [("series_", "series"), ("term_", "term"), ("snapshot_", "snapshot"),
                  ("elsa_", None)]

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
    nml_read cannot parse. Values are kept as written (quotes included, the
    comment removed); scalar() and listval() read them.
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
        cur[key] = strip_comment(val).strip().rstrip(",").strip()
    return groups, errors


def strip_comment(val):
    """The value without its inline comment (a ! outside quotes)."""
    q = None
    for i, c in enumerate(val):
        if q:
            if c == q:
                q = None
        elif c in "'\"":
            q = c
        elif c == "!":
            return val[:i]
    return val


def scalar(val):
    """A scalar value without its quotes."""
    return val.strip().strip("'\"")


def listval(val):
    """A list value ("a" "b", or unquoted words) as a list, "" entries dropped."""
    items = re.findall(r'"([^"]*)"|\'([^\']*)\'|([^\s,"\']+)', val)
    return [a or b or c for a, b, c in items if (a or b or c)]


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
            name = scalar(groups[ctrl].get(key, defaults["yelmo"].get(key, def_group)))
            out.setdefault(name, def_group)
    return out


def domain_suffixes(groups):
    """The domain suffixes of a file: "" and <sfx> of every &domain<sfx>."""
    return sorted({""} | {g[len("domain"):] for g in groups if g.startswith("domain")})


def yelmox_groups(groups, ydefaults):
    """User group -> template group of yelmox_defaults.nml for the forcing
    groups of a file ([forcing<sfx>], the series and terms they list, the
    snapshot backend [snapshots<sfx>], the transforms), and the named groups
    looked up through the one var_defs list ([forcing] var_defs) -> "forcing".
    The snapshot groups follow once the terms are resolved
    (snapshot_groups)."""
    out, stale, owner = {}, [], {}
    if "obm_forcing" in groups:
        out["obm_forcing"] = "obm_forcing"
    for g in groups:
        tmpl = library_template(g)
        if tmpl:
            out[g] = tmpl
    forcing = groups.get("forcing", {})
    for n in listval(forcing.get("series", ydefaults["forcing"].get("series", ""))):
        out.setdefault("series_" + n, "series")
        owner["series_" + n] = "forcing"
    for sfx in domain_suffixes(groups):
        fg = "forcing" + sfx
        if fg in groups:
            out[fg] = "forcing"
        for base in ("taper", "bias", "orbit", "ice_bnd"):
            if base + sfx in groups:
                out[base + sfx] = base
        for n in listval(groups.get(fg, {}).get("terms", ydefaults["forcing"].get("terms", ""))):
            out.setdefault("term_" + n, "term")
            owner["term_" + n] = "forcing"
        climate = scalar(groups.get("comps" + sfx, {}).get("climate", ""))
        if climate in ("snapshots", "rembo"):
            if "snapshots" + sfx in groups:
                out["snapshots" + sfx] = "snapshots"
            ref = scalar(groups.get(fg, {}).get("reference", ydefaults["forcing"]["reference"]))
            out.setdefault("snapshot_" + ref, "snapshot")
            owner["snapshot_" + ref] = "forcing"
    return out, stale, owner


def snapshot_groups(groups, xgroups, owner, resolved):
    """The snapshot groups the resolved snapshot terms name (snap0, snap1):
    added to xgroups/owner; True when any was new."""
    new = False
    for name, (src, params) in list(resolved.items()):
        if xgroups.get(name) != "term":
            continue
        if scalar(params.get("pattern", "snapshots")) != "snapshots":
            continue
        for key in ("snap0", "snap1"):
            sn = scalar(params.get(key, ""))
            if sn and ("snapshot_" + sn) not in xgroups:
                xgroups["snapshot_" + sn] = "snapshot"
                owner["snapshot_" + sn] = "forcing"
                new = True
    return new


def var_defs_files(path, listed, label, errors):
    """The par file itself (searched implicitly first), then the var_defs
    files of a group in order; a listed file that is missing is an error."""
    files = []
    for f in listed:
        # Relative to the run folder: the repo links (input, ...).
        cand = [os.path.join(REPO, f), os.path.join(os.path.dirname(path), f)]
        found = next((c for c in cand if os.path.isfile(c)), None)
        if found is None:
            errors.append("&%s: var_defs file not found: %s" % (label, f))
        else:
            files.append(found)
    return [path] + files


def library_groups(path, cache):
    """Group name -> keys (first copies) of a var_defs file, parsed once."""
    if path not in cache:
        cache[path] = first_copies(parse_par(path)[0])[0]
    return cache[path]


def resolve_named(path, groups, ygroups, owner, vd_files, cache):
    """Where each forcing group of a file is read from, as group -> (file,
    keys): the par file for its own groups; a group looked up through a
    var_defs list (owner) from the par file, else the first var_defs file
    that holds it. Also returns the errors for looked-up groups found in
    none."""
    out, errors = {}, []
    for name in ygroups:
        if name in groups:
            out[name] = (path, groups[name])
            continue
        if name not in owner:
            continue
        files = vd_files[owner[name]]
        found = next((f for f in files[1:] if name in library_groups(f, cache)), None)
        if found is None:
            errors.append("&%s (from &%s): group in neither the par file nor the var_defs files"
                          % (name, owner[name]))
        else:
            out[name] = (found, library_groups(found, cache)[name])
    return out, errors


def var_defs_errors(path, groups, ygroups, ydefaults, vd_files, resolved, cache):
    """Bindings of the snapshots (and [bias] dtas/fpr, [taper] weight, [ice_bnd]
    H_ice) to varslice groups that neither the par file nor the var_defs files
    hold. Also returns the bound groups found in the par file itself (read by
    varslice, not dead)."""
    errors, local = [], set()
    files = vd_files["forcing"]

    def check(label, refs):
        db = set()
        for f in files[1:]:
            db |= set(library_groups(f, cache))
        for key, ref in refs:
            if ref in groups:
                local.add(ref)
            elif ref not in db:
                errors.append("&%s %s: group %s in neither the par file nor the var_defs files"
                              % (label, key, ref))

    for name, tmpl in ygroups.items():
        if tmpl == "snapshot":
            if name not in resolved:
                continue
            keys = resolved[name][1]
            check(name, [(key, ref) for key in ("tas", "pr", "zs", "to", "so")
                         for ref in listval(keys.get(key, ""))])
        elif tmpl in ("bias", "taper", "ice_bnd"):
            g = groups[name]
            if scalar(g.get("active", ydefaults[tmpl]["active"])).lower() not in (
                    "true", ".true.", "t"):
                continue
            keys = {"bias": ("dtas", "fpr"), "taper": ("weight",), "ice_bnd": ("H_ice",)}[tmpl]
            check(name, [(key, scalar(g.get(key, ydefaults[tmpl].get(key, ""))))
                         for key in keys
                         if scalar(g.get(key, ydefaults[tmpl].get(key, "")))])
    return errors, local


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


def library_template(name):
    """Template group of a library group (by prefix), None if not checked,
    False if the prefix is not one of a library."""
    return next((t for p, t in LIBRARY_PREFIX if name.startswith(p)), False)


def check_library(groups, errors, ydefaults):
    """Return (errors, dead) for a library file: each group against its
    template by prefix."""
    for name, params in groups.items():
        tmpl = library_template(name)
        if not tmpl:
            continue
        unknown = [k for k in params if k not in ydefaults[tmpl]]
        if unknown:
            errors.append("&%s: unknown key(s) (not in &%s of yelmox defaults): %s"
                          % (name, tmpl, ", ".join(unknown)))
    return errors, []


def check_file(path, defaults, ydefaults, lits):
    """Return (errors, dead) for one par file or library file, or None if it
    is neither."""
    raw, errors = parse_par(path)
    groups, dup = first_copies(raw)
    if not control_blocks(groups):
        if groups and all(library_template(g) is not False for g in groups):
            return check_library(groups, errors + dup, ydefaults)
        return None
    errors += dup

    ygroups = yelmo_groups(groups, defaults)
    for name, def_group in ygroups.items():
        unknown = [k for k in groups.get(name, {}) if k not in defaults[def_group]]
        if unknown:
            errors.append("&%s: unknown key(s) (not in &%s of defaults): %s"
                          % (name, def_group, ", ".join(unknown)))

    xgroups, stale, owner = yelmox_groups(groups, ydefaults)
    # The one var_defs list, [forcing] var_defs (the par file first, then the
    # listed files); the snapshot groups follow the resolved terms.
    cache = {}
    vd_files = {"forcing": var_defs_files(path, listval(groups.get("forcing", {}).get(
        "var_defs", ydefaults["forcing"]["var_defs"])), "forcing", errors)}
    resolved, rs_errors = resolve_named(path, groups, xgroups, owner, vd_files, cache)
    if snapshot_groups(groups, xgroups, owner, resolved):
        resolved, rs_errors = resolve_named(path, groups, xgroups, owner, vd_files, cache)
    errors += rs_errors
    for name, def_group in xgroups.items():
        if name not in resolved:
            continue
        src, params = resolved[name]
        unknown = [k for k in params if k not in ydefaults[def_group]]
        if unknown:
            where = "" if src == path else " (in %s)" % os.path.relpath(src, REPO)
            errors.append("&%s%s: unknown key(s) (not in &%s of yelmox defaults): %s"
                          % (name, where, def_group, ", ".join(unknown)))
    vd_errors, local = var_defs_errors(path, groups, xgroups, ydefaults, vd_files,
                                       resolved, cache)
    errors += vd_errors

    dead = []
    if lits is not None:
        dead += ["&%s: group not read by the snapshot backend" % g for g in stale]
        for name, params in groups.items():
            if name in ygroups or name in xgroups or name in stale or name in local:
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
    ap.add_argument("--yelmox-defaults", default=YX_DEFAULTS,
                    help="yelmox forcing defaults file (default: input/yelmox_defaults.nml)")
    ap.add_argument("--dead", action="store_true",
                    help="also report keys of the other groups no source reads (errors)")
    args = ap.parse_args(argv[1:])

    loaded = []
    for path in (args.defaults, args.yelmox_defaults):
        if not os.path.isfile(path):
            print("ERROR: defaults file not found: %s" % path)
            return 2
        d, def_errors = parse_par(path)
        d, dup = first_copies(d)
        if def_errors or dup:
            print("ERROR: defaults file %s: %s" % (path, "; ".join(def_errors + dup)))
            return 2
        loaded.append(d)
    defaults, ydefaults = loaded

    lits = None
    if args.dead:
        missing = [d for d in SRC_DIRS if not os.path.isdir(os.path.join(REPO, d))]
        if missing:
            print("ERROR: --dead needs the source trees: %s" % ", ".join(missing))
            return 2
        lits = source_literals()
    files = args.files or sorted(f for g in PAR_GLOBS
                                 for f in glob.glob(os.path.join(REPO, g)))

    print("Defaults: %s, %s" % (os.path.relpath(args.defaults, REPO),
                                os.path.relpath(args.yelmox_defaults, REPO)))
    print()

    n_fail = 0
    for path in files:
        rel = os.path.relpath(path, REPO)
        result = check_file(path, defaults, ydefaults, lits)
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
