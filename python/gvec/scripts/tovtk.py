# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""tovtk.py - Save a VTK file with some default values from a GVEC output file."""

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
    prog="pygvec-tovtk", description="Output a VTK file for reading into paraview."
)

parser.add_argument("--rundir", type=Path, help="GVEC run directory", default=Path("."))

parser.add_argument(
    "--outputfile", "-o", type=Path, help="VTK output filename", default=Path("./out")
)

parser.add_argument(
    "--verbose",
    "-v",
    action="count",
    default=0,
    help="verbosity level: -v for info, -vv for debug",
)

parser.add_argument(
    "--ntheta",
    type=int,
    default=50,
    help="theta resolution",
)

parser.add_argument(
    "--nzeta",
    type=int,
    default=33,
    help="zeta resolution",
)

parser.add_argument(
    "--period",
    choices=["full", "single", "half"],
    default="full",
    help="zeta resolution",
)


def gvec_to_vtk(
    state: State,
    outputfile: Path,
    ntheta: int = 50,
    nzeta: int = 33,
    period: Literal["full", "single", "half"] = "full",
    verbose: bool = False,
):
    if period == "full":
        zeta_adjust = 1.0
    elif period == "single":
        zeta_adjust = state.nfp
    elif period == "half":
        zeta_adjust = state.nfp / 2.0

    theta = np.linspace(0, 2 * np.pi, ntheta)
    zeta = np.linspace(0, 2 * np.pi / zeta_adjust, nzeta)

    evaluations = state.evaluate(
        "X1",
        "X2",
        "LA",
        "pos",
        "B",
        "J",
        "mod_B",
        "p",
        "iota",
        "Phi",
        "chi",
        "I_tor",
        "I_pol",
        rho=20,
        theta=theta,
        zeta=zeta,
    )

    ev2vtk(outputfile, evaluations, not verbose)


def main(args: Sequence[str] | argparse.Namespace | None = None):
    if isinstance(args, argparse.Namespace):
        pass
    else:
        args = parser.parse_args(args)

    state = gvec.find_state(args.rundir)

    if args.verbose > 0:
        verbose = True
    else:
        verbose = False

    gvec_to_vtk(state, args.outputfile, args.ntheta, args.nzeta, args.period, verbose)


if __name__ == "__main":
    main()
