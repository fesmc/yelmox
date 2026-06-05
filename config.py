"""
Generate the top-level Makefile with the correct build configuration.

Usage:
    python config.py <machine> <compiler>     # e.g. macbook gfortran
    python config.py --file <legacy_config>    # e.g. --file macbook_gfortran

In the two-axis form, two fragments are assembled into the Makefile
template (config/Makefile), in this order:

    config/compilers/<compiler>.mk   compiler flags (FC, FFLAGS, DFLAGS, ...)
    config/machines/<machine>.mk     machine paths (NetCDF) and any overrides

The shared dependency wiring (config/common.mk) is `include`d by the
Makefile template itself, after the placeholder -- it is no longer
concatenated here. A machine fragment may still override a compiler
default (e.g. DFLAGS_NODEBUG), and common.mk references variables both
fragments set.

The --file form reproduces the legacy behaviour, dropping a single
per-host config file (kept in config/legacy/) into the template.
"""

import argparse
import os
import re
import sys
from glob import glob

CONFIG_DIR = "config"
COMPILERS_DIR = os.path.join(CONFIG_DIR, "compilers")
MACHINES_DIR = os.path.join(CONFIG_DIR, "machines")
LEGACY_DIR = os.path.join(CONFIG_DIR, "legacy")
TEMPLATE = os.path.join(CONFIG_DIR, "Makefile")
PLACEHOLDER = "<COMPILER_CONFIGURATION>"


def parse_makefile(filename="Makefile"):
    """Print the list of available `make` targets in a Makefile."""
    try:
        with open(filename, "r") as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return

    print("Available `make` commands (targets):")
    for line in lines:
        line = line.split("#")[0].strip()
        match = re.match(r"^([a-zA-Z0-9_-]+)\s*:\s*(.*)", line)
        if match and not line.startswith("\t"):
            print(f"  - {match.group(1)}")
    print("")


def available(directory):
    """Return sorted base names of the *.mk fragments in a directory."""
    return sorted(
        os.path.splitext(os.path.basename(p))[0]
        for p in glob(os.path.join(directory, "*.mk"))
    )


def read(path, kind):
    """Read a fragment, exiting with a helpful message if it is missing."""
    if not os.path.isfile(path):
        print(f"Error: {kind} config not found: {path}\n")
        if kind == "compiler":
            print("Available compilers:", ", ".join(available(COMPILERS_DIR)))
        elif kind == "machine":
            print("Available machines: ", ", ".join(available(MACHINES_DIR)))
        sys.exit(1)
    return open(path).read()


def main():
    parser = argparse.ArgumentParser(
        description="Generate the top-level Makefile build configuration.",
    )
    parser.add_argument("machine", nargs="?", help="Machine name, e.g. macbook")
    parser.add_argument("compiler", nargs="?", help="Compiler name, e.g. gfortran")
    parser.add_argument(
        "--file",
        metavar="LEGACY_CONFIG",
        help="Use a single legacy config file from config/legacy/ instead "
        "of assembling machine + compiler fragments.",
    )
    args = parser.parse_args()

    if args.file:
        config_path = os.path.join(LEGACY_DIR, args.file)
        if not os.path.isfile(config_path):
            print(f"Error: legacy config not found: {config_path}\n")
            print("Available legacy configs:", ", ".join(available(LEGACY_DIR)))
            sys.exit(1)
        compile_info = open(config_path).read()
        source_desc = config_path
    else:
        if not (args.machine and args.compiler):
            parser.error(
                "give both <machine> and <compiler> (or use --file). "
                "Try: python config.py macbook gfortran"
            )
        compiler_mk = read(
            os.path.join(COMPILERS_DIR, args.compiler + ".mk"), "compiler"
        )
        machine_mk = read(
            os.path.join(MACHINES_DIR, args.machine + ".mk"), "machine"
        )
        # common.mk is `include`d directly by the Makefile template; the
        # machine fragment may still override compiler defaults.
        compile_info = "\n".join([compiler_mk, machine_mk])
        source_desc = f"machine={args.machine}, compiler={args.compiler}"

    template = open(TEMPLATE).read()
    makefile = template.replace(PLACEHOLDER, compile_info)
    open("Makefile", "w").write(makefile)

    print(f"\nMakefile configuration complete for: {source_desc}\n")
    parse_makefile("Makefile")


if __name__ == "__main__":
    main()
