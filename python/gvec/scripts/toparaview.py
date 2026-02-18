# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""toparaview.py - Save a paraview VTS file with some default values from a GVEC output file."""

import argparse
import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Literal

import numpy as np

import gvec
from gvec import State
from gvec.util import logging_setup
from gvec.vtk import ev2vtk

parser = argparse.ArgumentParser(
    prog="pygvec-toparaview",
    description="Output a VTS file for reading into paraview.",
    formatter_class=argparse.RawTextHelpFormatter,
)

parser.add_argument("--rundir", type=Path, help="GVEC run directory", default=Path("."))

parser.add_argument(
    "--outputfile", "-o", type=Path, help="VTK output filename", default=Path("./out")
)

parser.add_argument(
    "--variables",
    type=str,
    default="X1, X2, LA, B, J, grad_rho, dp_dr",
    help="Variables to output to paraview file, as comma separated string. Default: 'X1, X2, LA, B, J, grad_rho, dp_dr'",
)

parser.add_argument(
    "--verbose",
    "-v",
    action="store_true",
    help="verbosity level: -v for debug",
)

parser.add_argument(
    "--nrho",
    type=int,
    default=11,
    help="rho resolution. Default: 11",
)

parser.add_argument(
    "--ntheta",
    type=int,
    default=41,
    help="theta resolution. Default: 41",
)

parser.add_argument(
    "--nzeta",
    type=int,
    default=51,
    help="zeta resolution per field period. Default: 51",
)

parser.add_argument(
    "--period",
    choices=["full", "single", "half"],
    default="single",
    help="Save data for a full, single or half period. Default: single",
)


def gvec_to_paraview(
    state: State,
    outputfile: Path,
    nrho: int = 11,
    ntheta: int = 41,
    nzeta: int = 51,
    period: Literal["full", "single", "half"] = "full",
    variables: list[str] = ["X1", "X2", "LA", "B", "J", "grad_rho", "dp_dr"],
    verbose: bool = False,
):
    if period == "full":
        zeta_end = 2.0 * np.pi
        nzeta = nzeta * state.nfp
    elif period == "single":
        zeta_end = 2.0 * np.pi / state.nfp
    elif period == "half":
        zeta_end = np.pi / state.nfp
        nzeta = nzeta // 2

    theta = np.linspace(0, 2 * np.pi, ntheta)
    zeta = np.linspace(0, zeta_end, nzeta)

    output_variables = ["pos", *variables]

    evaluations = state.evaluate(
        *output_variables,
        rho=20,
        theta=theta,
        zeta=zeta,
    )
    evaluations = evaluations[output_variables]

    ev2vtk(outputfile, evaluations, not verbose)


def main(args: Sequence[str] | argparse.Namespace | None = None):
    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)

    state = gvec.find_state(args.rundir)

    output_variables = args.variables.replace(" ", "").split(",")

    gvec_to_paraview(
        state,
        args.outputfile,
        args.nrho,
        args.ntheta,
        args.nzeta,
        args.period,
        output_variables,
        args.verbose,
    )


if __name__ == "__main":
    main()
