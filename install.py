#!/usr/bin/env python3
"""
Interactive installer for YelmoX.

Each component repository has its own install_* function that performs the
steps specific to that repo: clone, pick + apply config, plus any repo-local
extras (e.g. internal symlinks, `make install`, etc.). To add a new dependency
later, write a new install_<repo>(state) function and call it from main().

The script also writes a `.install.sh` that contains every bash command it
ran, so the install can be reproduced or copied into the docs without re-
answering any prompts.

Compilation of the static libraries is left to the user — the script prints
next-step build commands at the end.

Run from the yelmox repo root:
    python3 install.py
    python3 install.py --https      # clone via HTTPS instead of SSH
"""

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_CONFIG = "macbook_gfortran"
DEFAULT_HPC = "dkrz_levante"
DEFAULT_ACCOUNT = "ba1442"


# ---------------------------------------------------------------- shared state

@dataclass
class State:
    yelmox_root: Path
    install_dir: Path
    protocol: str
    default_cfg: str
    include_rembo: bool
    hpc: str
    account: str
    ice_data_path: str
    iso_data_path: str
    repo_paths: dict = field(default_factory=dict)
    pending_configs: list = field(default_factory=list)
    data_pending: list = field(default_factory=list)
    bash_log: list = field(default_factory=list)
    bash_cwd: Path = None


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

def log_section(state, title):
    state.bash_log.append("")
    state.bash_log.append(f"# --- {title} ---")


def _log_cd(state, cwd):
    if cwd is None:
        return
    cwd = Path(cwd)
    if cwd != state.bash_cwd:
        state.bash_log.append(f"cd {shlex.quote(str(cwd))}")
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


def clone(state, name, org, dest):
    if dest.exists():
        print(f"  - {name}: already at {dest}, skipping clone")
        # Still record the clone command so .install.sh is reproducible from scratch.
        url = (f"git@github.com:{org}/{name}.git" if state.protocol == "ssh"
               else f"https://github.com/{org}/{name}.git")
        log_cmd(state, ["git", "clone", url], cwd=dest.parent)
        return
    url = (f"git@github.com:{org}/{name}.git" if state.protocol == "ssh"
           else f"https://github.com/{org}/{name}.git")
    print(f"  $ git clone {url}   (in {dest.parent})")
    log_cmd(state, ["git", "clone", url], cwd=dest.parent)
    subprocess.run(["git", "clone", url], cwd=dest.parent, check=True)


def make_link(state, target, link_path):
    """Relative symlink link_path -> target. Skip if link_path already exists.
    Always logs the equivalent ln -s into the bash log so .install.sh is
    reproducible from scratch."""
    target = Path(target)
    link_path = Path(link_path)
    cwd = link_path.parent
    rel = os.path.relpath(target, cwd)
    log_cmd(state, ["ln", "-s", rel, link_path.name], cwd=cwd)

    if link_path.is_symlink() or link_path.exists():
        print(f"  - link {link_path.name} already exists in {cwd}, skipping")
        return
    try:
        link_path.symlink_to(rel)
        print(f"  + {link_path} -> {rel}")
    except ValueError:
        link_path.symlink_to(target)
        print(f"  + {link_path} -> {target}")


# ---------------------------------------------------------------- composed steps

def clone_into(state, name, org, dirname):
    """Clone <org>/<name> into <install_dir>/<dirname>, record in state, return path."""
    dest = state.install_dir / dirname
    clone(state, name, org, dest)
    state.repo_paths[name] = dest
    return dest


def do_config(state, repo_name, repo_dir):
    """Prompt for a config name (numbered list); run config.py or queue as pending."""
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
    """First yelmox step: copy .runme/runme_config to .runme_config and set hpc/account."""
    log_section(state, "yelmox runme_config")
    src = state.yelmox_root / ".runme" / "runme_config"
    dst = state.yelmox_root / ".runme_config"
    if not src.exists():
        print(f"  ! {src} missing — skipping runme_config setup")
        return
    if dst.exists():
        print(f"  - {dst} already exists, leaving as-is (delete it to regenerate)")
    else:
        print(f"  $ cp {src} {dst}")
        dst.write_bytes(src.read_bytes())
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

    # Equivalent sed commands in .install.sh (BSD/GNU-portable).
    log_raw(state,
            f"sed -i.bak -E 's/(\"hpc\"[[:space:]]*:[[:space:]]*\")[^\"]*(\")/"
            f"\\1{state.hpc}\\2/' .runme_config",
            cwd=state.yelmox_root)
    log_raw(state,
            f"sed -i.bak -E 's/(\"account\"[[:space:]]*:[[:space:]]*\")[^\"]*(\")/"
            f"\\1{state.account}\\2/' .runme_config",
            cwd=state.yelmox_root)
    log_raw(state, "rm .runme_config.bak", cwd=state.yelmox_root)

    print(f"  Note: edit {dst.name} by hand to change jobname, omp, mem, email, etc.")


def install_fesm_utils(state):
    log_section(state, "fesm-utils")
    dest = clone_into(state, "fesm-utils", "fesmc", "fesm-utils")
    do_config(state, "fesm-utils", dest)


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
    """yelmox itself: config + root-level links to siblings + external data links."""
    log_section(state, "yelmox (config + root links)")
    do_config(state, "yelmox", state.yelmox_root)

    # Root-level links to component repos
    make_link(state, state.repo_paths["fesm-utils"],   state.yelmox_root / "fesm-utils")
    make_link(state, state.repo_paths["yelmo"],        state.yelmox_root / "yelmo")
    make_link(state, state.repo_paths["FastIsostasy"], state.yelmox_root / "FastIsostasy")
    if state.include_rembo:
        make_link(state, state.repo_paths["rembo1"], state.yelmox_root / "rembo1")

    # External data dir links
    for label, path in [("ice_data", state.ice_data_path),
                        ("isostasy_data", state.iso_data_path)]:
        if path:
            p = Path(path).expanduser().resolve()
            if not p.exists():
                print(f"  ! {label} path {p} does not exist — creating link anyway")
            make_link(state, p, state.yelmox_root / label)
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


def collect_initial_inputs(yelmox_root, https_flag):
    print(f"YelmoX root: {yelmox_root}")
    print()
    print("This script will:")
    print("  1. Set up .runme_config (copy from .runme/runme_config, set hpc + account)")
    print("  2. Clone the YelmoX component repositories")
    print("  3. Run config.py in each with a machine config you pick")
    print("  4. Create the symlinks needed for the build")
    print("  5. Write the equivalent bash to .install.sh for reproducibility")
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

    print("\nHPCs defined in .runme/queues.json:")
    hpcs = load_hpc_names(yelmox_root)
    if hpcs:
        hpc = pick_numbered("Default HPC for .runme_config", hpcs, DEFAULT_HPC)
    else:
        hpc = ask("Default HPC for .runme_config", default=DEFAULT_HPC)
    account = ask("Account for .runme_config", default=DEFAULT_ACCOUNT)

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
        hpc=hpc,
        account=account,
        ice_data_path=ice_data_path,
        iso_data_path=iso_data_path,
    )


def write_install_sh(state):
    path = state.yelmox_root / ".install.sh"
    header = [
        "#!/usr/bin/env bash",
        "# Generated by install.py.",
        "# Reproduces a YelmoX install with the choices made in the interactive run.",
        "# Edit if your situation differs (e.g. machine config name, account, paths).",
        "set -euo pipefail",
    ]
    body = state.bash_log
    path.write_text("\n".join(header + body) + "\n")
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

    if state.data_pending:
        print(f"\nPending data links (run from {state.yelmox_root}):")
        for label in state.data_pending:
            print(f"  ln -s /path/to/{label} {label}")

    print(f"\n.runme_config: hpc={state.hpc}, account={state.account} set.")
    print(f"  Edit by hand to change jobname, omp, mem, email, etc.")

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
    print(f"  {step}. Compile yelmox:")
    print(f"       cd {state.yelmox_root}")
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

    print(f"\n=== Setting up yelmox runme_config ===")
    install_runme(state)

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

    print(f"\n=== Writing .install.sh ===")
    write_install_sh(state)

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
