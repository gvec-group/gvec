# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.core.state.util - utilities for working with the gvec.State class"""

# === Imports === #

import logging
from pathlib import Path

from ._class import State

# === Globals === #

logger = logging.getLogger("gvec.state")

# === Functions === #


def load_state(parameterfile: Path | str, statefile: Path | str):
    """Load a State object from a given parameterfile (.ini) and statefile (.dat)."""
    return State(parameterfile, statefile)


def find_state(rundir: Path | str | None = None):
    """Load a State object from a given run-directory. Use the latest statefile which is found."""
    if rundir is None:
        rundir = Path.cwd()
    else:
        rundir = Path(rundir)
    parameterfiles = list(rundir.glob("parameter*.ini"))
    if len(parameterfiles) > 1:
        # If there is more than one parameter file we check if a "final" set is present
        parameterfiles_final = [file for file in parameterfiles if "final" in file.name]
        if len(parameterfiles_final) == 1:
            logger.info("found more than one candidate parameterfile, with one '*final*' file")
            parameterfiles = parameterfiles_final
        else:
            raise ValueError(
                f"found more than one candidate parameterfile: {[file.name for file in parameterfiles]}"
            )
    elif len(parameterfiles) == 0:
        raise ValueError("no parameterfile found")
    logger.info(f"found parameterfile '{parameterfiles[0].name}'")

    statefiles = sorted(rundir.glob("*State*.dat"))
    projectnames = set([f.name.split("_State")[0] for f in statefiles])
    if len(statefiles) == 0:
        raise ValueError("no statefile found")
    if len(projectnames) > 1:
        raise ValueError(
            f"found statefiles for different projects: {projectnames}; cannot determine which one to load"
        )
    logger.info(f"found statefile '{statefiles[-1].name}'")

    return State(parameterfiles[0], statefiles[-1])


def find_states(rundir: Path | str | None = None):
    """LEGACY: Load a Sequence of State objects from a given run-directory."""
    if rundir is None:
        rundir = Path.cwd()
    else:
        rundir = Path(rundir)
    parameterfiles = list(rundir.glob("parameter*.ini"))
    if len(parameterfiles) > 1:
        raise ValueError(
            f"Found more than one candidate parameterfile: {[file.name for file in parameterfiles]}"
        )
    elif len(parameterfiles) == 0:
        raise ValueError("No parameterfile found.")
    statefiles = sorted(rundir.glob("*State*.dat"))
    return [State(parameterfiles[0], statefile) for statefile in statefiles]


def load_intermediates(rundir: Path | str | None = None):
    """Load a Sequence of State objects from a given run-directory."""
    if rundir is None:
        rundir = Path.cwd()
    else:
        rundir = Path(rundir)
    stages = sorted(rundir.glob("stage-*"))
    states_out = {}
    for stage in stages:
        stage_states = []
        parameterfiles = sorted(stage.glob("*_parameters_*.ini"))
        statefiles = sorted(stage.glob("*State*"))
        if len(parameterfiles) == 0:
            raise ValueError("No parameterfile found.")
        if len(statefiles) == 0:
            raise ValueError("No statefile found")
        for i_pf, params in enumerate(parameterfiles[:-1]):
            iteration_start = _get_file_iteration(params)
            iteration_end = _get_file_iteration(parameterfiles[i_pf + 1])
            stage_states += _associate_statefiles(
                statefiles, params, iteration_start, iteration_end
            )
        iteration_start = _get_file_iteration(parameterfiles[-1])
        stage_states += _associate_statefiles(
            statefiles, parameterfiles[-1], iteration_start, int(1e7)
        )
        states_out[stage.name] = stage_states
    return states_out


def _associate_statefiles(
    statefiles: list[Path], params: Path, iteration_start: int, iteration_end: int
):
    """Helper function for find_states. Create a list of States associated with a set of parameters

    Parameters
    ----------
    statefiles : list[Path]
        List of possible statefiles
    params : Path
        Parameter file to associate to
    iteration_start : int
        minimum number of iterations
    iteration_end : int
        maximum number of iterations

    Returns
    -------
    list[State]
        List of states associated with the given parameter file.
    """
    states_out = []
    for sf in statefiles:
        iteration_sf = int(sf.name.split("_")[-1].split(".")[0])
        if iteration_sf >= iteration_start and iteration_sf < iteration_end:
            states_out.append(State(params, sf))
        elif iteration_sf < iteration_start:
            continue
        else:
            break
    return states_out


def _get_file_iteration(file_in: Path):
    """Get the iteration number from a gvec state or parameter file.

    Parameters
    ----------
    file_in : Path
        Filepath to the state or parameter file

    Returns
    -------
    int
        Iteration number
    """
    # TODO: use some fancy re for this?
    return int(file_in.name.split("_")[-1].split(".")[0])
