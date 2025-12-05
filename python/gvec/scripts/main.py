# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""The pyGVEC executable"""

# === Imports === #

import platform
from pathlib import Path
import logging
import argparse
from collections.abc import Sequence

import gvec
from gvec.scripts import cas3d, convert, gist, quasr, run

# === Arguments === #

parser = argparse.ArgumentParser(
    prog="pygvec",
    description=f"GVEC: a flexible 3D MHD equilibrium solver\npyGVEC v{gvec.__version__}",
)
subparsers = parser.add_subparsers(
    title="mode",
    description="which mode/subcommand to run",
    dest="mode",
    metavar="MODE",
)
parser.add_argument(
    "-V",
    "--version",
    action="version",
    version=gvec.util.version_info(),
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

convert_parser = subparsers.add_parser(
    "convert-params",
    help="convert the GVEC parameterfile between different formats",
    formatter_class=convert.parser.formatter_class,
    description=convert.parser.description,
    parents=[convert.parser],
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

cas3d_parser = subparsers.add_parser(
    "to-cas3d",
    help="convert a GVEC state to a CAS3D compatible input file",
    description=cas3d.parser.description,
    parents=[cas3d.parser],
    add_help=False,
)

gist_parser = subparsers.add_parser(
    "to-gist",
    help="convert a GVEC state to a GENE-GIST compatible input file",
    description=gist.parser.description,
    parents=[gist.parser],
    add_help=False,
)


# === Script === #


def main(args: Sequence[str] | argparse.Namespace | None = None):
    gvec.util.logging_setup()

    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)

    # --- run GVEC scripts --- #
    if args.mode == "run":
        return run.main(args)

    elif args.mode == "convert-params":
        return convert.main(args)

    elif args.mode == "to-cas3d":
        return cas3d.main(args)

    elif args.mode == "to-gist":
        return gist.main(args)

    elif args.mode == "load-quasr":
        return quasr.main(args)


if __name__ == "__main__":
    exit(main())
