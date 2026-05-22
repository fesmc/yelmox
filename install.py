#!/usr/bin/env python3
"""
Interactive installer for YelmoX.

Each component repository has its own install_* function that performs the
steps specific to that repo: clone, configure for the current system, plus
any repo-local extras (e.g. internal symlinks, `make install`, etc.). To add
a new dependency later, write a new install_<repo>(state) function and call
it from main().

The script also writes a `.install.sh` that contains every bash command it
ran, so the install can be reproduced or copied into the docs without re-
answering any prompts.

Compilation of the static libraries is left to the user — the script prints
next-step build commands at the end.

Run from the yelmox repo root:
    python3 install.py
    python3 install.py -d clone-https        # clone via HTTPS instead of SSH
    python3 install.py -d no                 # repos already under yelmox; just link
    python3 install.py --no-config           # skip every config step (.runme + per-repo)
    python3 install.py --config-default mymachine   # change the default machine config
    python3 install.py --overwrite           # re-clone repos (existing -> outdated-repos/)
"""

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_CONFIG = "macbook_gfortran"
DEFAULT_HPC = "dkrz_levante"
DEFAULT_ACCOUNT = "ba1442"

# Mapping from clone-name to bash variable name used in .install.sh.
REPO_VAR_NAMES = {
    "fesm-utils":   "fesm_utils",
    "yelmo":        "yelmo",
    "FastIsostasy": "fastiso",
    "coordinates":  "coord",
    "rembo1":       "rembo1",
}


# ---------------------------------------------------------------- shared state

@dataclass
class State:
    yelmox_root: Path
    install_dir: Path                # ignored when download == "no"
    download: str                    # "clone-ssh" | "clone-https" | "no"
    do_config_steps: bool            # False when --no-config is passed
    default_cfg: str
    include_rembo: bool
    hpc: str
    account: str
    ice_data_path: str
    iso_data_path: str
    overwrite: bool = False
    repo_paths: dict = field(default_factory=dict)
    pending_configs: list = field(default_factory=list)
    data_pending: list = field(default_factory=list)
    bash_log: list = field(default_factory=list)
    bash_cwd: Path = None
    # Ordered maps of bash variables emitted at the top of .install.sh.
    # path_vars: name -> absolute Path (defined in clone order, used to substitute
    #   paths in cd/ln commands).
    # string_vars: name -> string (e.g. hpc, account).
    path_vars: dict = field(default_factory=dict)
    string_vars: dict = field(default_factory=dict)


# ---------------------------------------------------------------- io helpers

def ask(prompt, default=None, allow_empty=False):
    suffix = ""
    if default:
        suffix = f" [{default}]"
    elif allow_empty:
        suffix = " [skip]"
    while True:
        try:
            ans = input(f"{prompt}{suffix}: ").strip()
        except EOFError:
            ans = ""
        if ans:
            return ans
        if default is not None:
            return default
        if allow_empty:
            return ""
        print("  (please enter a value)")


def ask_yn(prompt, default=False):
    d = "Y/n" if default else "y/N"
    while True:
        try:
            ans = input(f"{prompt} [{d}]: ").strip().lower()
        except EOFError:
            ans = ""
        if not ans:
            return default
        if ans in ("y", "yes"):
            return True
        if ans in ("n", "no"):
            return False
        print("  (answer y or n)")


def pick_numbered(prompt, options, default):
    """Show a numbered list; accept either a number or a literal name.
    Returns the chosen string. A name not in `options` is returned as-is
    (caller can treat that as a custom value)."""
    print(f"  Options:")
    for i, opt in enumerate(options, start=1):
        marker = " *" if opt == default else ""
        print(f"    {i:2d}) {opt}{marker}")
    suggested = default if default in options else (options[0] if options else None)
    raw = ask(f"  {prompt} (number, name, or a new name for custom)",
              default=suggested)
    if raw.isdigit():
        idx = int(raw)
        if 1 <= idx <= len(options):
            return options[idx - 1]
        print(f"  ! number {idx} out of range; treating as literal name")
    return raw


# ---------------------------------------------------------------- bash log

def path_to_bash(state, p):
    """Render an absolute Path as a bash expression, using a defined path
    variable when possible (longest exact match wins; otherwise longest prefix).
    Falls back to a quoted absolute path."""
    p = Path(p)
    if not p.is_absolute():
        return shlex.quote(str(p))
    # Exact match wins
    for name, pp in state.path_vars.items():
        if pp == p:
            return f'"${name}"'
    # Otherwise longest-prefix wins
    items = sorted(state.path_vars.items(),
                   key=lambda kv: len(str(kv[1])), reverse=True)
    for name, pp in items:
        try:
            rel = p.relative_to(pp)
        except ValueError:
            continue
        return f'"${name}/{rel}"'
    return shlex.quote(str(p))


def _bash_rhs_for(state, current_name, abs_path):
    """RHS for `current_name=<...>` in the .install.sh header. Considers only
    variables registered before current_name so dependencies resolve left-to-right."""
    abs_path = Path(abs_path)
    for name, pp in state.path_vars.items():
        if name == current_name:
            break
        if pp == abs_path:
            return f'"${name}"'
        try:
            rel = abs_path.relative_to(pp)
        except ValueError:
            continue
        return f'"${name}/{rel}"'
    return shlex.quote(str(abs_path))


def register_path_var(state, name, p):
    """Register a path under a bash variable name (idempotent)."""
    state.path_vars[name] = Path(p)


def log_section(state, title):
    state.bash_log.append("")
    state.bash_log.append(f"# --- {title} ---")


def _log_cd(state, cwd):
    if cwd is None:
        return
    cwd = Path(cwd)
    if cwd != state.bash_cwd:
        state.bash_log.append(f"cd {path_to_bash(state, cwd)}")
        state.bash_cwd = cwd


def log_cmd(state, cmd, cwd=None):
    _log_cd(state, cwd)
    state.bash_log.append(" ".join(shlex.quote(str(c)) for c in cmd))


def log_raw(state, line, cwd=None):
    _log_cd(state, cwd)
    state.bash_log.append(line)


# ---------------------------------------------------------------- primitives

def run(state, cmd, cwd=None, check=True):
    where = f"   (in {cwd})" if cwd else ""
    print(f"  $ {' '.join(str(c) for c in cmd)}{where}")
    log_cmd(state, cmd, cwd=cwd)
    return subprocess.run(cmd, cwd=cwd, check=check)


def list_configs(repo_dir):
    cfg = repo_dir / "config"
    if not cfg.is_dir():
        return []
    out = []
    for p in sorted(cfg.iterdir()):
        if not p.is_file():
            continue
        name = p.name
        if name.startswith(".") or name.startswith("Makefile") or name.endswith(".mk"):
            continue
        out.append(name)
    return out


def _clone_url(state, name, org):
    if state.download == "clone-https":
        return f"https://github.com/{org}/{name}.git"
    return f"git@github.com:{org}/{name}.git"


def clone(state, name, org, dest):
    """Clone <org>/<name> into <dest>. When download == 'no', do not clone:
    only verify that <dest> exists (as a real dir or symlink) and record a
    comment in .install.sh."""
    if state.download == "no":
        if not (dest.exists() or dest.is_symlink()):
            raise FileNotFoundError(
                f"{name}: -d no requires {dest} to already exist (as dir or symlink)")
        print(f"  - {name}: using existing {dest} (no clone)")
        log_raw(state,
                f"# {name}: assumed already present at {path_to_bash(state, dest)} "
                f"(provide via git clone, symlink, or extracted archive)")
        return

    url = _clone_url(state, name, org)

    # --overwrite: move any existing copy aside before cloning fresh. Done
    # silently (not logged to .install.sh) — it's a one-shot cleanup step, not
    # part of the reproducible bash sequence.
    if state.overwrite and (dest.exists() or dest.is_symlink()):
        outdated = state.install_dir / "outdated-repos"
        moved = outdated / dest.name
        if moved.exists() or moved.is_symlink():
            raise FileExistsError(
                f"{moved} already exists from a previous --overwrite run; "
                f"remove it manually before re-running install.py --overwrite.")
        outdated.mkdir(exist_ok=True)
        print(f"  - {name}: moving existing {dest} -> {moved}")
        shutil.move(str(dest), str(moved))

    log_cmd(state, ["git", "clone", url], cwd=dest.parent)
    if dest.exists():
        print(f"  - {name}: already at {dest}, skipping clone")
        return
    print(f"  $ git clone {url}   (in {dest.parent})")
    subprocess.run(["git", "clone", url], cwd=dest.parent, check=True)


def make_link(state, target, link_path, absolute=False):
    """Symlink link_path -> target. Skip if link_path already exists.
    By default the on-disk link is relative (portable); pass absolute=True to
    force an absolute link regardless of how target was provided. In .install.sh
    we emit the target via a bash variable when one matches, since variables make
    the generated script easier to edit by hand."""
    target = Path(target)
    link_path = Path(link_path)
    cwd = link_path.parent
    rel = os.path.relpath(target, cwd)

    _log_cd(state, cwd)
    state.bash_log.append(
        f"ln -s {path_to_bash(state, target)} {shlex.quote(link_path.name)}")

    if link_path.is_symlink() or link_path.exists():
        print(f"  - {link_path.name} already exists in {cwd}, skipping")
        return
    if absolute:
        link_path.symlink_to(target)
        print(f"  + {link_path} -> {target}")
        return
    try:
        link_path.symlink_to(rel)
        print(f"  + {link_path} -> {rel}")
    except ValueError:
        link_path.symlink_to(target)
        print(f"  + {link_path} -> {target}")


# ---------------------------------------------------------------- composed steps

def clone_into(state, name, org, dirname):
    """Resolve the path for <name>, clone if needed, record in state, return path.
    When download == 'no', the canonical path is yelmox_root/<dirname> (which may
    be a real dir or a symlink the user already arranged); otherwise it is
    install_dir/<dirname>."""
    if state.download == "no":
        dest = state.yelmox_root / dirname
    else:
        dest = state.install_dir / dirname
    # Register the repo's path under a bash variable before logging the clone,
    # so any subsequent log line (including the -d no `# assumed at ...` comment)
    # can refer to it as $name.
    if name in REPO_VAR_NAMES:
        register_path_var(state, REPO_VAR_NAMES[name], dest)
    state.repo_paths[name] = dest
    clone(state, name, org, dest)
    return dest


def do_config(state, repo_name, repo_dir):
    """Prompt for a config name (numbered list); run config.py or queue as pending.
    No-op when --no-config was passed."""
    if not state.do_config_steps:
        return
    configs = list_configs(repo_dir)
    print(f"\n  --- {repo_name} ({repo_dir}) ---")
    if not configs:
        print(f"  (no config/ directory found, skipping)")
        return
    cfg = pick_numbered("Choose config", configs, state.default_cfg)
    cfg_path = repo_dir / "config" / cfg

    if cfg in configs:
        run(state, ["python3", "config.py", f"config/{cfg}"], cwd=repo_dir)
        return

    # Custom config: copy a template if file doesn't yet exist, then queue as pending.
    if not cfg_path.exists():
        tmpl = None
        if state.default_cfg in configs:
            tmpl = repo_dir / "config" / state.default_cfg
        elif configs:
            tmpl = repo_dir / "config" / configs[0]
        if tmpl is not None:
            print(f"  {repo_name}: creating {cfg_path.name} from template {tmpl.name}")
            cfg_path.write_bytes(tmpl.read_bytes())
            log_raw(state,
                    f"cp {shlex.quote(f'config/{tmpl.name}')} {shlex.quote(f'config/{cfg}')}",
                    cwd=repo_dir)
        else:
            print(f"  ! {repo_name}: no template available; create {cfg_path} yourself")
    else:
        print(f"  {repo_name}: {cfg_path.name} already exists, leaving as-is")
    state.pending_configs.append((repo_name, repo_dir, cfg, cfg_path))
    print(f"  (pending) edit {cfg_path}, then run config.py")
    log_raw(state,
            f"# TODO: edit config/{cfg} for your machine, then:",
            cwd=repo_dir)
    log_raw(state, f"# python3 config.py config/{cfg}", cwd=repo_dir)


# ---------------------------------------------------------------- per-repo install funcs
#
# To add a new dependency repo:
#   def install_<name>(state):
#       dest = clone_into(state, "<gh-repo>", "<gh-org>", "<dirname>")
#       do_config(state, "<gh-repo>", dest)
#       # any repo-specific extras here, e.g.:
#       #   run(state, ["make", "install"], cwd=dest)
#       #   make_link(state, state.repo_paths["fesm-utils"], dest / "fesm-utils")
# Then add a call to install_<name>(state) from main() in the right order.

def install_runme(state):
    """Create the local .runme_config via `runme --config`, then set hpc/account.
    `runme --config` copies .runme/runme_config to .runme_config (the canonical
    first-clone step). No-op when --no-config was passed."""
    if not state.do_config_steps:
        return
    log_section(state, "yelmox runme_config")
    src = state.yelmox_root / ".runme" / "runme_config"
    dst = state.yelmox_root / ".runme_config"
    if not src.exists():
        print(f"  ! {src} missing — skipping runme_config setup")
        return
    if dst.exists():
        print(f"  - {dst} already exists, leaving as-is (delete it and run `runme --config` to regenerate)")
        log_raw(state, "# .runme_config already present; to regenerate: runme --config", cwd=state.yelmox_root)
    elif shutil.which("runme"):
        # On a fresh clone there is no local .runme_config yet, so `runme --config`
        # copies the template silently (no interactive overwrite prompt).
        run(state, ["runme", "--config"], cwd=state.yelmox_root)
    else:
        print(f"  ! `runme` not on PATH — copying {src.name} to {dst.name} directly")
        print(f"  $ cp {src} {dst}")
        dst.write_bytes(src.read_bytes())
        log_raw(state, "# runme not on PATH; equivalent of `runme --config`:", cwd=state.yelmox_root)
        log_raw(state, "cp .runme/runme_config .runme_config", cwd=state.yelmox_root)

    # Patch hpc and account in place, preserving formatting.
    text = dst.read_text()
    new_text = re.sub(r'("hpc"\s*:\s*")[^"]*(")',
                      rf'\g<1>{state.hpc}\g<2>', text)
    new_text = re.sub(r'("account"\s*:\s*")[^"]*(")',
                      rf'\g<1>{state.account}\g<2>', new_text)
    if new_text != text:
        dst.write_text(new_text)
        print(f"  set hpc={state.hpc}, account={state.account} in {dst.name}")

    # Register hpc and account as string variables; the sed commands below
    # reference them by name so they're easy to edit at the top of .install.sh.
    state.string_vars["hpc"] = state.hpc
    state.string_vars["account"] = state.account

    # Equivalent sed commands in .install.sh (BSD/GNU-portable). Use bash
    # double-quotes around the sed expression so ${hpc}/${account} expand.
    log_raw(state,
            'sed -i.bak -E "s/(\\"hpc\\"[[:space:]]*:[[:space:]]*\\")[^\\"]*(\\")/'
            '\\1${hpc}\\2/" .runme_config',
            cwd=state.yelmox_root)
    log_raw(state,
            'sed -i.bak -E "s/(\\"account\\"[[:space:]]*:[[:space:]]*\\")[^\\"]*(\\")/'
            '\\1${account}\\2/" .runme_config',
            cwd=state.yelmox_root)
    log_raw(state, "rm .runme_config.bak", cwd=state.yelmox_root)

    print(f"  Note: edit {dst.name} by hand to change jobname, omp, mem, email, etc.")


def install_runme_pkg(state):
    """Install fesmc/runme via pip (provides the `runme` command used to stage,
    run, and submit simulations and ensembles). Skip if `runme` is already on
    PATH. No-op when --no-config was passed."""
    if not state.do_config_steps:
        return
    log_section(state, "runme (runme command)")

    url = "git+https://github.com/fesmc/runme"

    if shutil.which("runme"):
        print(f"  - `runme` already available at {shutil.which('runme')}, skipping install")
        log_raw(state, "# runme already installed: `runme` command found on PATH")
        log_raw(state, f"# To reinstall: pip install {url}")
        return

    print("  `runme` command not found on PATH.")
    if not ask_yn("  Install runme via pip from fesmc/runme", default=True):
        print("  Skipping runme install — `runme` command will not be available.")
        log_raw(state, "# TODO: install runme manually:")
        log_raw(state, f"# pip install {url}")
        return

    run(state, [sys.executable, "-m", "pip", "install", url])

    runme_path = shutil.which("runme")
    if runme_path:
        print(f"  + `runme` is now available at {runme_path}")
    else:
        print("  ! pip install succeeded but `runme` is not on PATH.")
        print("    The Python bin dir is probably missing from PATH.")
        print("    Add to ~/.bashrc or ~/.profile:")
        print('      PATH=${PATH}:${HOME}/.local/bin')
        print('      export PATH')
        log_raw(state, "# `runme` not on PATH after install — add Python user bin to PATH:")
        log_raw(state, '# PATH=${PATH}:${HOME}/.local/bin')
        log_raw(state, "# export PATH")


def install_fesm_utils(state):
    log_section(state, "fesm-utils")
    dest = clone_into(state, "fesm-utils", "fesmc", "fesm-utils")
    # fesm-utils is special: config.py and config/ live in the utils/ subdir,
    # not at the repo root. The generated Makefile also lives there. The
    # yelmox-side $(FESMUTILSROOT) points at fesm-utils/utils for this reason.
    do_config(state, "fesm-utils", dest / "utils")
    # Reminder in .install.sh about the one-time external-deps build.
    # Left as commented examples because the choice is machine/compiler-specific.
    log_raw(state, "# One-time external deps (LIS + FFTW). Slow (~10-30 min).")
    log_raw(state, "# Pick the script for your machine + compiler, e.g.:")
    log_raw(state, '# (cd "$fesm_utils" && ./install.sh)                  # generic')
    log_raw(state, '# (cd "$fesm_utils" && ./install_dkrz.sh ifx)         # DKRZ Levante')
    log_raw(state, '# (cd "$fesm_utils" && ./install_awi.sh  gfortran)    # AWI albedo')
    log_raw(state, '# (cd "$fesm_utils" && ./install_pik.sh  ifx)         # PIK HPC')


def install_coordinates(state):
    log_section(state, "coordinates")
    dest = clone_into(state, "coordinates", "fesmc", "coordinates")
    do_config(state, "coordinates", dest)


def install_yelmo(state):
    log_section(state, "yelmo")
    dest = clone_into(state, "yelmo", "fesmc", "yelmo")
    do_config(state, "yelmo", dest)
    make_link(state, state.repo_paths["fesm-utils"], dest / "fesm-utils")


def install_fastisostasy(state):
    log_section(state, "FastIsostasy")
    dest = clone_into(state, "FastIsostasy", "palma-ice", "FastIsostasy")
    do_config(state, "FastIsostasy", dest)
    make_link(state, state.repo_paths["fesm-utils"], dest / "fesm-utils")


def install_rembo1(state):
    log_section(state, "rembo1")
    dest = clone_into(state, "rembo1", "alex-robinson", "rembo1")
    do_config(state, "rembo1", dest)
    libs = dest / "libs"
    libs.mkdir(exist_ok=True)
    make_link(state, state.repo_paths["coordinates"], libs / "coordinates")


def install_yelmox(state):
    """yelmox itself: config + root-level links to siblings. External data
    links (ice_data, isostasy_data) are handled separately by setup_data_links."""
    log_section(state, "yelmox (config + root links)")
    do_config(state, "yelmox", state.yelmox_root)

    # Root-level links to component repos. With -d no the repo paths *are* the
    # yelmox-root entries, so the links are no-ops; skip them.
    if state.download != "no":
        make_link(state, state.repo_paths["fesm-utils"],   state.yelmox_root / "fesm-utils")
        make_link(state, state.repo_paths["yelmo"],        state.yelmox_root / "yelmo")
        make_link(state, state.repo_paths["FastIsostasy"], state.yelmox_root / "FastIsostasy")
        if state.include_rembo:
            # coordinates also needs a yelmox-root link: the parent Makefile's
            # COORDROOT = coordinates points here for the coord-static build.
            make_link(state, state.repo_paths["coordinates"], state.yelmox_root / "coordinates")
            make_link(state, state.repo_paths["rembo1"],      state.yelmox_root / "rembo1")


def setup_data_links(state):
    """Prompt for ice_data and isostasy_data paths and create symlinks under
    yelmox/. Either may be skipped (the summary will remind the user to wire
    them up later). Run after the per-repo installs and before runme_config —
    these are local-filesystem paths, separate from HPC settings."""
    log_section(state, "external data links")
    state.ice_data_path = ask("Path to ice_data (Enter to skip and link later)",
                              default="", allow_empty=True)
    state.iso_data_path = ask("Path to isostasy_data (Enter to skip and link later)",
                              default="", allow_empty=True)

    if state.ice_data_path:
        register_path_var(state, "ice_data_src",
                          Path(state.ice_data_path).expanduser().resolve())
    if state.iso_data_path:
        register_path_var(state, "isostasy_data_src",
                          Path(state.iso_data_path).expanduser().resolve())

    for label, path in [("ice_data", state.ice_data_path),
                        ("isostasy_data", state.iso_data_path)]:
        if path:
            p = Path(path).expanduser().resolve()
            if not p.exists():
                print(f"  ! {label} path {p} does not exist — creating link anyway")
            make_link(state, p, state.yelmox_root / label, absolute=True)
        else:
            state.data_pending.append(label)


# ---------------------------------------------------------------- orchestration

def load_hpc_names(yelmox_root):
    qf = yelmox_root / ".runme" / "queues.json"
    if not qf.exists():
        return []
    try:
        data = json.loads(qf.read_text())
    except json.JSONDecodeError as e:
        print(f"  ! could not parse {qf}: {e}")
        return []
    return sorted(data.keys())


def collect_initial_inputs(yelmox_root, download, do_config_steps,
                           default_cfg, overwrite):
    print(f"YelmoX root: {yelmox_root}")
    print(f"Download mode: {download}    Config steps: {'on' if do_config_steps else 'OFF (--no-config)'}")
    if overwrite:
        print("Overwrite mode: ON (existing repos move to install_dir/outdated-repos/)")
    print()
    print("This script will:")
    if download == "no":
        print("  1. Use the component repos already present under yelmox/")
    else:
        print(f"  1. Clone the YelmoX component repositories ({download})")
    if do_config_steps:
        print(f"  2. Configure each repo for the current system (default: {default_cfg})")
    print("  3. Create the symlinks needed for the build")
    if do_config_steps:
        print("  4. Install runme (fesmc/runme) via pip if `runme` is not already on PATH")
        print("  5. Set up .runme_config (via `runme --config`, set hpc + account)")
    print("  6. Write the equivalent bash to .install.sh for reproducibility")
    print()

    if download == "no":
        install_dir = yelmox_root  # not really used; repos live under yelmox/
    else:
        install_dir = Path(ask("Directory to clone component repos into",
                               default=str(Path.cwd()))).expanduser().resolve()
        install_dir.mkdir(parents=True, exist_ok=True)

    include_rembo = ask_yn("Include REMBO support (clones coordinates + rembo1)",
                           default=False)

    state = State(
        yelmox_root=yelmox_root,
        install_dir=install_dir,
        download=download,
        do_config_steps=do_config_steps,
        default_cfg=default_cfg,
        include_rembo=include_rembo,
        hpc=DEFAULT_HPC,
        account=DEFAULT_ACCOUNT,
        ice_data_path="",
        iso_data_path="",
        overwrite=overwrite,
    )

    # Register the bash variables that will appear at the top of .install.sh.
    # Order matters: later vars can be defined in terms of earlier ones. The
    # ice_data_src / isostasy_data_src vars are registered later from inside
    # setup_data_links, after the user has provided the paths.
    register_path_var(state, "yelmox_root", yelmox_root)
    if download != "no":
        register_path_var(state, "install_dir", install_dir)

    return state


def collect_runme_inputs(state):
    """Prompt for runme_config hpc/account at the end of the install, after
    clones and per-repo configs are done. Asks first whether the defaults are
    acceptable, so users on a local machine can skip the HPC-specific prompts
    entirely. No-op when --no-config was passed."""
    if not state.do_config_steps:
        return
    print()
    print(f"runme hpc configuration defaults ok? (hpc={DEFAULT_HPC}, account={DEFAULT_ACCOUNT})")
    if ask_yn("  Use defaults", default=True):
        return
    print("\nHPCs defined in .runme/queues.json:")
    hpcs = load_hpc_names(state.yelmox_root)
    if hpcs:
        state.hpc = pick_numbered("Default HPC for .runme_config", hpcs, DEFAULT_HPC)
    else:
        state.hpc = ask("Default HPC for .runme_config", default=DEFAULT_HPC)
    state.account = ask("Account for .runme_config", default=DEFAULT_ACCOUNT)


def cleanup_outdated_repos(state):
    """If --overwrite moved existing repos into install_dir/outdated-repos,
    ask at the very end whether to delete that folder. Not logged to
    .install.sh — this is a one-shot cleanup step, not part of the
    reproducible bash sequence."""
    if not state.overwrite:
        return
    outdated = state.install_dir / "outdated-repos"
    if not outdated.exists():
        return
    print()
    if ask_yn(f"Delete {outdated} (pre-overwrite repos)", default=False):
        shutil.rmtree(outdated)
        print(f"  removed {outdated}")
    else:
        print(f"  kept {outdated}")


def write_install_sh(state):
    path = state.yelmox_root / ".install.sh"
    out = [
        "#!/usr/bin/env bash",
        "# Generated by install.py.",
        "# Reproduces a YelmoX install with the choices made in the interactive run.",
        "# Edit the variables below to retarget paths or change hpc/account.",
        "set -euo pipefail",
    ]
    if state.path_vars:
        out.append("")
        out.append("# === Paths ===")
        for name, p in state.path_vars.items():
            out.append(f"{name}={_bash_rhs_for(state, name, p)}")
    if state.string_vars:
        out.append("")
        out.append("# === Settings ===")
        for name, v in state.string_vars.items():
            out.append(f"{name}={shlex.quote(v)}")
    out.extend(state.bash_log)
    path.write_text("\n".join(out) + "\n")
    try:
        path.chmod(0o755)
    except OSError:
        pass
    print(f"  + wrote {path}")


def print_summary(state):
    print()
    print("=" * 60)
    print("Install complete.")
    print("=" * 60)

    print(f"\nWrote .install.sh — equivalent bash sequence for the choices made.")
    print(f"  {state.yelmox_root / '.install.sh'}")

    if state.pending_configs:
        print("\nPending custom configs — edit each file, then run config.py:")
        for repo_name, repo_dir, cfg, cfg_path in state.pending_configs:
            print(f"  - {repo_name}:")
            print(f"      edit {cfg_path}")
            print(f"      (cd {repo_dir} && python3 config.py config/{cfg})")

    if state.do_config_steps:
        print(f"\n.runme_config, current settings:")
        print(f"    hpc     = {state.hpc}")
        print(f"    account = {state.account}")
        print(f"  (edit by hand to change them, as well as jobname, omp, mem, email, etc.)")
    else:
        print(f"\nConfig steps skipped (--no-config). .runme_config and per-repo")
        print(f"  config.py invocations were NOT run.")

    fesm = state.repo_paths["fesm-utils"]
    print("\nNext steps:")
    print(f"  1. Prep fesm-utils external deps (LIS + FFTW) — ONE-TIME, slow (~10-30 min):")
    print(f"       cd {fesm}")
    print(f"       # Pick the script for your machine + compiler, e.g.:")
    print(f"       # ./install.sh                  # generic")
    print(f"       # ./install_dkrz.sh ifx         # DKRZ Levante")
    print(f"       # ./install_awi.sh  gfortran    # AWI albedo")
    print(f"       # ./install_pik.sh  ifx         # PIK HPC")
    print(f"     If you skip this, `make yelmox` will abort early with a")
    print(f"     LIS/FFTW-not-built error (no half-built state).")
    extras = "fesm-utils, yelmo, FastIsostasy"
    if state.include_rembo:
        extras += ", coordinates, rembo1"
    print(f"  2. Compile yelmox (builds {extras} too):")
    print(f"       cd {state.yelmox_root}")
    print(f"       make clean")
    print(f"       make yelmox          # default build")
    if state.include_rembo:
        print(f"       make yelmox_rembo    # REMBO-coupled build")
    if state.data_pending:
        print(f"  3. Set up data links — needed at runtime:")
        for label in state.data_pending:
            print(f"       ln -s /path/to/{label} ./")
        print(f"")


def main():
    parser = argparse.ArgumentParser(description="Interactive YelmoX installer.")
    parser.add_argument("-d", "--download",
                        choices=["clone-ssh", "clone-https", "no", "none"],
                        default="clone-ssh",
                        help="How to obtain component repos. "
                             "'clone-ssh' (default) and 'clone-https' do a git clone. "
                             "'no'/'none' assumes each repo is already available at "
                             "yelmox/<repo> (as a dir or symlink) and only sets up "
                             "internal links + configs.")
    parser.add_argument("--no-config",
                        action="store_true",
                        help="Skip all configuration steps: do not copy/edit "
                             ".runme_config and do not run python3 config.py in "
                             "any repo. Useful when you want to re-set up links "
                             "without touching existing configuration.")
    parser.add_argument("--config-default",
                        default=DEFAULT_CONFIG,
                        metavar="NAME",
                        help=f"Default machine config name used as the per-repo "
                             f"suggestion (default: {DEFAULT_CONFIG}). No "
                             f"interactive prompt — set it here to skip that "
                             f"question entirely.")
    parser.add_argument("--overwrite",
                        action="store_true",
                        help="Re-clone every component repo. Any existing copy "
                             "is moved to install_dir/outdated-repos/<name> "
                             "before the fresh clone, and you are asked at the "
                             "end whether to delete that folder. Not compatible "
                             "with '-d no'/'-d none'. The mv/rm steps are NOT "
                             "recorded in .install.sh.")
    args = parser.parse_args()
    download = "no" if args.download == "none" else args.download

    if args.overwrite and download == "no":
        parser.error("--overwrite is not compatible with -d no/-d none "
                     "(no clone happens in that mode).")

    yelmox_root = Path(__file__).resolve().parent
    state = collect_initial_inputs(yelmox_root, download,
                                   do_config_steps=not args.no_config,
                                   default_cfg=args.config_default,
                                   overwrite=args.overwrite)

    if state.download == "no":
        print(f"\n=== Using component repos already under {state.yelmox_root} ===")
    else:
        print(f"\n=== Installing components into {state.install_dir} ===")
    install_fesm_utils(state)
    if state.include_rembo:
        install_coordinates(state)
    install_yelmo(state)
    install_fastisostasy(state)
    if state.include_rembo:
        install_rembo1(state)

    print(f"\n=== Configuring yelmox ===")
    install_yelmox(state)

    print(f"\n=== Installing runme (`runme` command) ===")
    install_runme_pkg(state)

    print(f"\n=== External data links ===")
    setup_data_links(state)

    if state.do_config_steps:
        print(f"\n=== Setting up yelmox runme_config ===")
        collect_runme_inputs(state)
        install_runme(state)

    print(f"\n=== Writing .install.sh ===")
    write_install_sh(state)

    print_summary(state)
    cleanup_outdated_repos(state)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nAborted.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"\nCommand failed: {' '.join(str(c) for c in e.cmd)}")
        sys.exit(e.returncode or 1)
