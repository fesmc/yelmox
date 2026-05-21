#!/usr/bin/env python3
"""
Interactive installer for YelmoX.

Each component repository has its own install_* function that performs the
steps specific to that repo: clone, pick + apply config, plus any repo-local
extras (e.g. internal symlinks, `make install`, etc.). To add a new dependency
later, write a new install_<repo>(state) function and call it from main().

Compilation of the static libraries is left to the user — the script prints
next-step build commands at the end.

Run from the yelmox repo root:
    python3 install.py
    python3 install.py --https      # clone via HTTPS instead of SSH
"""

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_CONFIG = "macbook_gfortran"


# ---------------------------------------------------------------- shared state

@dataclass
class State:
    yelmox_root: Path
    install_dir: Path
    protocol: str
    default_cfg: str
    include_rembo: bool
    ice_data_path: str
    iso_data_path: str
    repo_paths: dict = field(default_factory=dict)
    pending_configs: list = field(default_factory=list)
    data_pending: list = field(default_factory=list)


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


def run(cmd, cwd=None, check=True):
    where = f"   (in {cwd})" if cwd else ""
    print(f"  $ {' '.join(str(c) for c in cmd)}{where}")
    return subprocess.run(cmd, cwd=cwd, check=check)


# ---------------------------------------------------------------- primitives

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


def clone(name, org, dest, protocol):
    if dest.exists():
        print(f"  - {name}: already at {dest}, skipping clone")
        return
    url = (f"git@github.com:{org}/{name}.git" if protocol == "ssh"
           else f"https://github.com/{org}/{name}.git")
    run(["git", "clone", url, str(dest)])


def make_link(target, link_path):
    """Relative symlink link_path -> target. Skip if link_path already exists."""
    if link_path.is_symlink() or link_path.exists():
        print(f"  - link {link_path.name} already exists in {link_path.parent}, skipping")
        return
    try:
        rel = os.path.relpath(target, link_path.parent)
        link_path.symlink_to(rel)
        print(f"  + {link_path} -> {rel}")
    except ValueError:
        link_path.symlink_to(target)
        print(f"  + {link_path} -> {target}")


# ---------------------------------------------------------------- composed steps

def clone_into(state, name, org, dirname):
    """Clone <org>/<name> into <install_dir>/<dirname>, record in state, return path."""
    dest = state.install_dir / dirname
    clone(name, org, dest, state.protocol)
    state.repo_paths[name] = dest
    return dest


def do_config(state, repo_name, repo_dir):
    """Prompt for a config name; run config.py or queue as a pending custom config."""
    configs = list_configs(repo_dir)
    print(f"\n  --- {repo_name} ({repo_dir}) ---")
    if not configs:
        print(f"  (no config/ directory found, skipping)")
        return
    print(f"  Available configs:")
    for c in configs:
        marker = " *" if c == state.default_cfg else "  "
        print(f"  {marker} {c}")
    suggested = state.default_cfg if state.default_cfg in configs else configs[0]
    cfg = ask("  Choose config (or type a new name for a custom one)", default=suggested)
    cfg_path = repo_dir / "config" / cfg

    if cfg in configs:
        run(["python3", "config.py", f"config/{cfg}"], cwd=repo_dir)
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
        else:
            print(f"  ! {repo_name}: no template available; create {cfg_path} yourself")
    else:
        print(f"  {repo_name}: {cfg_path.name} already exists, leaving as-is")
    state.pending_configs.append((repo_name, repo_dir, cfg, cfg_path))
    print(f"  (pending) edit {cfg_path}, then run config.py")


# ---------------------------------------------------------------- per-repo install funcs
#
# To add a new dependency repo:
#   def install_<name>(state):
#       dest = clone_into(state, "<gh-repo>", "<gh-org>", "<dirname>")
#       do_config(state, "<gh-repo>", dest)
#       # any repo-specific extras here, e.g.:
#       #   run(["make", "install"], cwd=dest)
#       #   make_link(state.repo_paths["fesm-utils"], dest / "fesm-utils")
# Then add a call to install_<name>(state) from main() in the right order.

def install_fesm_utils(state):
    dest = clone_into(state, "fesm-utils", "fesmc", "fesm-utils")
    do_config(state, "fesm-utils", dest)


def install_coordinates(state):
    dest = clone_into(state, "coordinates", "fesmc", "coordinates")
    do_config(state, "coordinates", dest)


def install_yelmo(state):
    dest = clone_into(state, "yelmo", "fesmc", "yelmo")
    do_config(state, "yelmo", dest)
    make_link(state.repo_paths["fesm-utils"], dest / "fesm-utils")


def install_fastisostasy(state):
    dest = clone_into(state, "FastIsostasy", "palma-ice", "FastIsostasy")
    do_config(state, "FastIsostasy", dest)
    make_link(state.repo_paths["fesm-utils"], dest / "fesm-utils")


def install_rembo1(state):
    dest = clone_into(state, "rembo1", "alex-robinson", "rembo1")
    do_config(state, "rembo1", dest)
    libs = dest / "libs"
    libs.mkdir(exist_ok=True)
    make_link(state.repo_paths["coordinates"], libs / "coordinates")


def install_yelmox(state):
    """yelmox itself: config + root-level links to siblings + external data links."""
    do_config(state, "yelmox", state.yelmox_root)

    # Root-level links to component repos
    make_link(state.repo_paths["fesm-utils"],   state.yelmox_root / "fesm-utils")
    make_link(state.repo_paths["yelmo"],        state.yelmox_root / "yelmo")
    make_link(state.repo_paths["FastIsostasy"], state.yelmox_root / "FastIsostasy")
    if state.include_rembo:
        make_link(state.repo_paths["rembo1"], state.yelmox_root / "rembo1")

    # External data dir links
    for label, path in [("ice_data", state.ice_data_path),
                        ("isostasy_data", state.iso_data_path)]:
        if path:
            p = Path(path).expanduser().resolve()
            if not p.exists():
                print(f"  ! {label} path {p} does not exist — creating link anyway")
            make_link(p, state.yelmox_root / label)
        else:
            state.data_pending.append(label)


# ---------------------------------------------------------------- orchestration

def collect_initial_inputs(yelmox_root, https_flag):
    print(f"YelmoX root: {yelmox_root}")
    print()
    print("This script will:")
    print("  1. Clone the YelmoX component repositories")
    print("  2. Run config.py in each with a machine config you pick")
    print("  3. Create the symlinks needed for the build")
    print("Compilation is left to you. Next-step build commands print at the end.")
    print()

    install_dir = Path(ask("Directory to clone component repos into",
                           default=str(yelmox_root.parent))).expanduser().resolve()
    install_dir.mkdir(parents=True, exist_ok=True)

    include_rembo = ask_yn("Include REMBO support (clones coordinates + rembo1)",
                           default=False)

    protocol = "ssh"
    if https_flag:
        protocol = "https"
    elif ask_yn("Use HTTPS for git clone (default is SSH)", default=False):
        protocol = "https"

    default_cfg = ask("Default machine config name (used as suggestion per repo)",
                      default=DEFAULT_CONFIG)

    ice_data_path = ask("Path to ice_data (Enter to skip and link later)",
                        default="", allow_empty=True)
    iso_data_path = ask("Path to isostasy_data (Enter to skip and link later)",
                        default="", allow_empty=True)

    return State(
        yelmox_root=yelmox_root,
        install_dir=install_dir,
        protocol=protocol,
        default_cfg=default_cfg,
        include_rembo=include_rembo,
        ice_data_path=ice_data_path,
        iso_data_path=iso_data_path,
    )


def print_summary(state):
    print()
    print("=" * 60)
    print("Install complete.")
    print("=" * 60)

    if state.pending_configs:
        print("\nPending custom configs — edit each file, then run config.py:")
        for repo_name, repo_dir, cfg, cfg_path in state.pending_configs:
            print(f"  - {repo_name}:")
            print(f"      edit {cfg_path}")
            print(f"      (cd {repo_dir} && python3 config.py config/{cfg})")

    if state.data_pending:
        print(f"\nPending data links (run from {state.yelmox_root}):")
        for label in state.data_pending:
            print(f"  ln -s /path/to/{label} {label}")

    fesm = state.repo_paths["fesm-utils"]
    print("\nNext steps (build):")
    step = 1
    print(f"  {step}. Build fesm-utils:")
    print(f"       cd {fesm}")
    print(f"       # ./install_<machine>.sh <compiler>   # see install_*.sh for options")
    print(f"       cd utils && make clean && make fesmutils-static")
    step += 1
    if state.include_rembo:
        coord = state.repo_paths['coordinates']
        print(f"  {step}. Build coordinates:")
        print(f"       cd {coord} && make clean && make coord-static")
        step += 1
    print(f"  {step}. Set runme config and compile yelmox:")
    print(f"       cd {state.yelmox_root}")
    print(f"       cp .runme/runme_config .runme_config   # then edit for your system")
    print(f"       make clean")
    print(f"       make yelmox          # default build")
    if state.include_rembo:
        print(f"       make yelmox_rembo    # REMBO-coupled build (uses rembo1 + coordinates)")


def main():
    parser = argparse.ArgumentParser(description="Interactive YelmoX installer.")
    parser.add_argument("--https", action="store_true",
                        help="Clone repos via HTTPS instead of SSH.")
    args = parser.parse_args()

    yelmox_root = Path(__file__).resolve().parent
    state = collect_initial_inputs(yelmox_root, args.https)

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

    print_summary(state)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nAborted.")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"\nCommand failed: {' '.join(str(c) for c in e.cmd)}")
        sys.exit(e.returncode or 1)
