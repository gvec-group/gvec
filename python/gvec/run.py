# ============================================================================================================================== #
# Copyright (C) 2024 Robert Babin <robert.babin@ipp.mpg.de>
# Copyright (C) 2024 Florian Hindenlang <hindenlang@gmail.com>

#
# This file is part of GVEC. GVEC is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.
#
# GVEC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
#
# You should have received a copy of the GNU General Public License along with GVEC. If not, see <http://www.gnu.org/licenses/>.
# ============================================================================================================================== #
"""pygvec run gvec from python"""

from . import _fgvec
from .lib import modgvec_py_run as _run
from .lib import modgvec_py_binding as _binding

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
