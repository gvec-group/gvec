# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""The pyGVEC executable"""

# === Imports === #

import platform
import time
from pathlib import Path
import shutil
import logging
import argparse
import re
from typing import Mapping, Sequence
from datetime import datetime
import copy

import numpy as np
import xarray as xr
import yaml
import tomlkit  # also supports writing (compared to tomllib)

import gvec
from gvec.scripts import to_cas3d, run

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
    version=gvec.__version_str__,
)

# --- convert parameterfile --- #

convert_parser = subparsers.add_parser(
    "convert-params",
    help="convert the GVEC parameterfile between different formats",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="Convert GVEC parameterfiles between different formats.\n"
    "The INI (classical) parameter files do not support stages or the current constraint!\nAlso the formatting is lost upon conversion.",
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
    description=to_cas3d.parser.description,
    parents=[to_cas3d.parser],
    add_help=False,
)

# === Script === #


def main(args: Sequence[str] | argparse.Namespace | None = None):
    logging.basicConfig(level=logging.WARNING)  # show warnings and above as normal
    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)

    # --- run GVEC --- #
    if args.mode == "run":
        run.main(args)

    # --- convert parameterfile --- #
    elif args.mode == "convert-params":
        parameters = gvec.util.read_parameters(args.input)
        gvec.util.write_parameters(parameters, args.output)

    # --- other scripts --- #
    elif args.mode == "to-cas3d":
        to_cas3d.main(args)


if __name__ == "__main__":
    main()
