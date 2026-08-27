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
from .fixed_boundary import FixedBoundaryRun
from .free_boundary import FreeBoundaryRun


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
    logger = logging.getLogger("gvec")
    if loglevel is not None:
        logger.setLevel(loglevel)
    else:
        logger.setLevel("WARNING")

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

    has_coilset = False
    # check if coil-set is specified in the parameters
    if "free_boundary" in params:
        if "coil_set" in params["free_boundary"]:
            has_coilset = True

    # make sure only one coil-set is specified
    if has_coilset and coil_set is not None:
        raise ValueError(
            "coil_set was handed as an input to run but is also present in the parameters! Either remove the coil_set from the specified parameters or do NOT hand it explicitly to gvec.run(...)."
        )
    elif has_coilset:
        coil_set = parameters["free_boundary"]["coil_set"]

    match coil_set:
        case CoilSet():
            do_free_boundary = True
            parameters["free_boundary"]["coil_set"] = "gvec-CoilSet"
        case str() | Path():
            do_free_boundary = True
            parameters["free_boundary"]["coil_set"] = str(coil_set)
            coil_set = CoilSet.load(coil_set)
        case None:
            do_free_boundary = False
        case _:
            if hasattr(coil_set, "eval_B") and callable(coil_set.eval_B):
                do_free_boundary = True
                parameters["free_boundary"]["coil_set"] = "custom field evaluation"
            else:
                raise TypeError(
                    f"Unknown coil_set type, expected either 'str', 'Path', 'CoilSet' or 'None' but got {type(coil_set)} which does not have a method 'eval_B(...)'"
                )

    # Run the case
    with gvec.util.chdir(runpath):
        if do_free_boundary:
            logger.info("Initializing free-boundary run")
            if keep_intermediates in ["stages", "all"]:
                save_fb_updates = True
            else:
                save_fb_updates = False
            run_instance = FreeBoundaryRun(
                params=params,
                coil_set=coil_set,
                state=restartstate,
                redirect_gvec_stdout=redirect_gvec_stdout,
                quiet=quiet,
                parameter_format=parameter_format,
                save_fb_updates=save_fb_updates,
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
