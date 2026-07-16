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
        final_found = False
        for param_file in parameterfiles:
            if "final" in param_file.name:
                parameterfiles = [parameterfiles[0]]
                final_found = True
        if not final_found:
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
    """Load a Sequence of State objects from a given run-directory."""
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
