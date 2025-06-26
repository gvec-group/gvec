# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""The pyGVEC executable"""

# === Imports === #

import platform
from pathlib import Path
import argparse
from collections.abc import Sequence

import gvec
from gvec.scripts import cas3d, run, quasr

# === Arguments === #

parser = argparse.ArgumentParser(
    prog="pygvec",
    description=f"GVEC: a flexible 3D MHD equilibrium solver\npyGVEC v{gvec.__version__}",
)
subparsers = parser.add_subparsers(
    title="mode",
    description="which mode/subcommand to run",
    dest="mode",
)
parser.add_argument(
    "-V",
    "--version",
    action="version",
    version=f"pyGVEC v{gvec.__version__} from {Path(gvec.__file__).parent} (python {platform.python_version()})",
)

# --- convert parameterfile --- #

convert_parser = subparsers.add_parser(
    "convert-params",
    help="convert the GVEC parameterfile between different formats",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Convert GVEC parameterfiles between different formats.\n"
    "The INI (classical) parameter files do not support stages or the current optimization!\nAlso the formatting is lost upon conversion.",
)
convert_parser.add_argument(
    "input",
    type=Path,
    help="input GVEC parameterfile",
)
convert_parser.add_argument(
    "output",
    type=Path,
    help="output GVEC parameterfile",
)

# --- scripts --- #

run_parser = subparsers.add_parser(
    "run",
    help="run GVEC (with stages)",
    formatter_class=run.parser.formatter_class,
    description=run.parser.description,
    parents=[run.parser],
    add_help=False,
)

cas3d_parser = subparsers.add_parser(
    "to-cas3d",
    help="convert a GVEC state to a CAS3D compatible input file",
    description=cas3d.parser.description,
    parents=[cas3d.parser],
    add_help=False,
)

quasr_parser = subparsers.add_parser(
    "load-quasr",
    help=quasr.parser.description,
    description=quasr.parser.description,
    parents=[quasr.parser],
    add_help=False,
    usage=quasr.parser.usage,
)

# === Script === #


def main(args: Sequence[str] | argparse.Namespace | None = None):
    gvec.util.logging_setup()

    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)

    # --- run GVEC --- #
    if args.mode == "run":
        return run.main(args)

    # --- convert parameterfile --- #
    elif args.mode == "convert-params":
        parameters = gvec.util.read_parameters(args.input)
        gvec.util.write_parameters(parameters, args.output)

    # --- other scripts --- #
    elif args.mode == "to-cas3d":
        return cas3d.main(args)

    elif args.mode == "load-quasr":
        return quasr.main(args)


if __name__ == "__main__":
    exit(main())
