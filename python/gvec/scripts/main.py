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
    help="input GVEC or VMEC parameterfile",
)
convert_parser.add_argument(
    "output",
    type=Path,
    nargs="?",
    help="output GVEC parameterfile",
    default="parameter.yaml",
)
convert_parser.add_argument(
    "--vmec",
    action="store_true",
    help="input parameterfile is a VMEC namelist",
)
convert_parser.add_argument(
    "-x",
    "--flip",
    choices=["t", "theta", "z", "zeta", "b", "both"],
    help="flip the coordinates in the specified direction(s), possible values are: t/theta, z/zeta, b/both",
    metavar="FLIP",
)
verbosity = convert_parser.add_mutually_exclusive_group()
verbosity.add_argument(
    "-v",
    "--verbose",
    action="count",
    default=0,
    help="verbosity level: -v for info, -vv for debug",
)
verbosity.add_argument("-q", "--quiet", action="store_true", help="suppress output")

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
    logger = logging.getLogger("gvec")

    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)

    # --- run GVEC --- #
    if args.mode == "run":
        return run.main(args)

    # --- convert parameterfile --- #
    elif args.mode == "convert-params":
        if args.quiet:
            logging.disable()
        elif args.verbose >= 2:
            logger.setLevel(logging.DEBUG)
        elif args.verbose == 1:
            logger.setLevel(logging.INFO)
        logger.debug(f"parsed args: {args}")

        if args.vmec:
            try:
                import f90nml
            except ImportError as e:
                logger.debug(f"Caught exception: {e}")
                logger.error("reading VMEC namelists requires 'f90nml' to be installed.")
            with open(args.input, "r") as file:
                content = file.read()
            content = content.strip()
            if content.endswith("&END"):
                content = content[:-4]
            nml = f90nml.reads(content)["indata"]
            parameters = gvec.util.parameters_from_vmec(nml, str(args.input))
        else:
            parameters = gvec.util.read_parameters(args.input)

        if args.flip is None:
            if not gvec.util.check_boundary_direction(parameters):
                logger.info("Input boundary is left-handed, flipping theta.")
                parameters = gvec.util.flip_parameters_theta(parameters)
        else:
            if args.flip[0] in "tb":
                parameters = gvec.util.flip_parameters_theta(parameters)
            if args.flip[0] in "zb":
                parameters = gvec.util.flip_parameters_zeta(parameters)

        if not gvec.util.check_boundary_direction(parameters):
            logger.warning("Output boundary is left-handed!")

        gvec.util.write_parameters(parameters, args.output)

    # --- other scripts --- #
    elif args.mode == "to-cas3d":
        return cas3d.main(args)

    elif args.mode == "load-quasr":
        return quasr.main(args)


if __name__ == "__main__":
    exit(main())
