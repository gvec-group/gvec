# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""run gvec from python"""

from . import _fgvec
from ._fgvec import modgvec_py_run as _run
from ._fgvec import modgvec_py_binding as _binding

from pathlib import Path


def run(
    parameterfile: str | Path,
    restartfile: str | Path | None = None,
    MPIcomm: int | None = None,
    stdout_path: str | Path | None = "stdout.txt",
):
    """
    Run gvec from python

    Parameters
    ----------
    parameterfile : str
        Path to / name of parameter file
    restartfile : str
        Path to / name of GVEC restart file, optional
    MPIcomm : int
        MPI communicator, optional (default in GVEC (if compiled with MPI) is MPI_COMM_WORLD)
    stdout_path : str
        Path to / name of file to redirect the standard output of GVEC. Optional, default is "stdout.txt".
        If set to None, stdout is not redirected
    """

    if stdout_path is not None:
        _binding.redirect_abort()
        _binding.redirect_stdout(str(stdout_path))

    if not Path(parameterfile).exists():
        raise FileNotFoundError(f"Parameter file {parameterfile} does not exist.")
    if restartfile is not None:
        if not Path(restartfile).exists():
            raise FileNotFoundError(f"Restart file {restartfile} does not exist.")

    _run.start_rungvec(str(parameterfile), restartfile_in=restartfile, comm_in=MPIcomm)
    # check for success
    if stdout_path is not None:
        with open(stdout_path) as stdout:
            if "GVEC SUCESSFULLY FINISHED!" not in stdout.read():
                raise RuntimeError(f"GVEC failed (see {stdout_path})")
