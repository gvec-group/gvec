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
from gvec.scripts import cas3d, gist as gist, run, quasr

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
    choices=["auto", "none", "tor", "pol", "both"],
    default="auto",
    help="flip the coordinates in the specified direction(s), possible values are: 'auto' (default), 'none', 'tor', 'pol', 'both'.",
    metavar="FLIP",
)
convert_parser.add_argument(
    "-t",
    "--shift-theta",
    choices=["auto", "pi", "0"],
    default="auto",
    help="shift the theta coordinate origin, possible values are 'auto' (default), 'pi' or '0'.",
    metavar="SHIFT",
)
convert_parser.add_argument(
    "--stellsym",
    action="store_true",
    help="enforce stellarator symmetry in the output parameterfile",
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

gist_parser = subparsers.add_parser(
    "to-gist",
    help="convert a GVEC state to a GENE-GIST compatible input file",
    description=gist.parser.description,
    parents=[gist.parser],
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
                logger.debug(f"caught exception: {e}")
                logger.error("reading VMEC namelists requires 'f90nml' to be installed.")
                return
            with open(args.input, "r") as file:
                content = file.read()
            content = content.strip()
            if content[-4:].lower() == "&end":
                content = content[:-4]
            nml = f90nml.reads(content)["indata"]
            parameters = gvec.util.parameters_from_vmec(nml, args.input.name)
            if args.flip == "auto":
                parameters = gvec.util.flip_parameters_zeta(parameters)
        else:
            parameters = gvec.util.read_parameters(args.input)

        if args.stellsym:
            parameters["X1_sin_cos"] = "_cos_"
            parameters["X2_sin_cos"] = "_sin_"
            parameters["LA_sin_cos"] = "_sin_"
            if "X1_b_sin" in parameters:
                del parameters["X1_b_sin"]
            if "X2_b_cos" in parameters:
                del parameters["X2_b_cos"]
            if "LA_b_cos" in parameters:
                del parameters["LA_b_cos"]
            if "X1_a_sin" in parameters:
                del parameters["X1_a_sin"]
            if "X2_a_cos" in parameters:
                del parameters["X2_a_cos"]

        if args.flip == "auto":
            if not gvec.util.check_boundary_direction(parameters):
                logger.info("input boundary is left-handed, flipping theta")
                parameters = gvec.util.flip_parameters_theta(parameters)
        elif args.flip in ["pol", "both"]:
            parameters = gvec.util.flip_parameters_theta(parameters)
        elif args.flip in ["tor", "both"]:
            parameters = gvec.util.flip_parameters_zeta(parameters)

        if not gvec.util.check_boundary_direction(parameters):
            logger.warning("output boundary is left-handed (increases clockwise)")

        if args.shift_theta == "auto":
            if "X1_sin_cos" not in parameters:
                logger.debug(
                    "cannot automatically shift theta origin: 'X1_sin_cos' not in parameters"
                )
            elif parameters["X1_sin_cos"] != "_cos_":
                logger.debug(
                    "cannot automatically shift theta origin: X1 is not stellarator symmetric"
                )
            elif "X1_b_cos" not in parameters:
                logger.debug(
                    "cannot automatically shift theta origin: 'X1_b_cos' not in parameters"
                )
            elif parameters["X1_b_cos"].get((1, 0)) < 0:
                logger.info(
                    "input boundary has theta origin at the 'inboard' side, shifting theta by pi"
                )
                parameters = gvec.util.shift_boundary_theta_pi(parameters)
        elif args.shift_theta == "pi":
            parameters = gvec.util.shift_boundary_theta_pi(parameters)

        gvec.util.write_parameters(parameters, args.output)

    # --- other scripts --- #
    elif args.mode == "to-cas3d":
        return cas3d.main(args)

    elif args.mode == "to-gist":
        return gist.main(args)

    elif args.mode == "load-quasr":
        return quasr.main(args)


if __name__ == "__main__":
    exit(main())
