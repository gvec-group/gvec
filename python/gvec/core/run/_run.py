# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""run gvec from python"""

import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Literal

import gvec
from gvec.core.state import State
from gvec.coils import CoilSet
from ._fixed_boundary_run import FixedBoundaryRun
from ._free_boundary_run import FreeBoundaryRun


def run(
    parameters: Mapping | Path | str,
    restartstate: Path | State | None = None,
    coil_set: CoilSet | Path | str | None = None,
    runpath: Path | str | None = None,
    redirect_gvec_stdout: bool = True,
    quiet: bool = False,
    parameter_format: Literal["toml", "yaml"] = "toml",
    keep_intermediates: Literal["all", "stages"] | None = None,
    loglevel: Literal["WARNING", "INFO", "DEBUG"] | None = None,
):
    """Run GVEC with the provided parameters.

    Parameters
    ----------
    parameters : Mapping | Path | str
        GVEC parameter dictionary or parameter file.
    restartstate : Path | State | None, optional
        Path to restart file or restart State object, by default None
    runpath : Path | str | None, optional
        Path to the directory where GVEC is executed. Overwrites existing directories or cerates a new directory if it does not exist.
        The default is None.
    redirect_gvec_stdout : bool, optional
        Whether to redirect GVEC's stdout. The default is True.
    quiet : bool, optional
        Whether to suppress all output. The default is False.
    parameter_format : Literal["toml", "yaml"], optional
        Format of the parameter file automatically generated if `picard_current="auto"`, by default "toml"
    keep_intermediates : Literal["all", "stages"] | None, optional
        Whether to keep intermediate results of GVEC. With `"all"`, all intermediate results are kept. With `"stages"`,
        only the final restarts from each stage are kept. With `None`, all intermediate results are deleted. The default is None.
    loglevel : Literal["WARNING", "INFO", "DEBUG"] | None, optional
        Set the loglevel for the gvec logger.

    Returns
    -------
    gvec.Run
        Run-object containing history and current state of the GVEC run.
    """

    # logger setup
    gvec.util.logging_setup()
    logger = logging.getLogger("gvec.run")
    if loglevel is not None:
        logger.setLevel(loglevel)

    # rundirectory setup
    if runpath is None:
        runpath = Path.cwd()
    else:
        runpath = Path(runpath)

        if not runpath.exists():
            runpath.mkdir()
            logger.info(f"created run directory {runpath}")
        else:
            logger.info(f"run directory {runpath} is overwritten")

    # transform parameters into hierarchical cidict.
    if isinstance(parameters, (Path, str)):
        parameters = Path(parameters)
        params = gvec.util.read_parameters(parameters)
    else:
        params = gvec.util.stack_parameters(parameters)

    match coil_set:
        case CoilSet():
            do_free_boundary = True
        case str() | Path():
            do_free_boundary = True
            coil_set = CoilSet.load(coil_set)
        case None:
            do_free_boundary = False
        case _:
            raise TypeError(
                f"Unknown coil_set type, expected either 'str', 'Path', 'CoilSet' or 'None' but got {type(coil_set)}"
            )

    # Run the case
    with gvec.util.chdir(runpath):
        if do_free_boundary:
            logger.info("Initializing free-boundary run")
            run_instance = FreeBoundaryRun(
                params=params,
                coil_set=coil_set,
                state=restartstate,
                redirect_gvec_stdout=redirect_gvec_stdout,
                quiet=quiet,
                parameter_format=parameter_format,
            )
        else:
            logger.info("Initializing fixed-boundary run")
            run_instance = FixedBoundaryRun(
                params=params,
                state=restartstate,
                redirect_gvec_stdout=redirect_gvec_stdout,
                quiet=quiet,
                parameter_format=parameter_format,
            )

        run_instance.run(keep_intermediates=keep_intermediates)

    return run_instance
