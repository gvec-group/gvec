# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""gvec.core.state._class - the gvec.State class"""

# === Imports === #

import logging
import tempfile
import warnings
from collections.abc import Mapping
from pathlib import Path

import xarray as xr

import gvec.lib
import gvec.util

from . import _binding
from ._binding import flush_stdout

# === Globals === #

logger = logging.getLogger("gvec.state")

# === State Class === #


class State:
    """A class for encapsulating the 'state' of the GVEC library with all relevant parameters and variables."""

    # === Constructor & Destructor === #

    def __init__(
        self,
        parameterfile: str | Path,
        statefile: str | Path | None = None,
    ):
        self.parameterfile: Path = Path(parameterfile).absolute()
        self.statefile: Path | None = Path(statefile).absolute() if statefile else None

        self._stdout: tempfile.NamedTemporaryFile | None = None
        self._children: list[gvec.lib.Modgvec_Sfl_Boozer.t_sfl_boozer] = []
        self.parameters: gvec.util.CaseInsensitiveDict = None

        if not self.parameterfile.exists():
            raise FileNotFoundError(f"Parameter file {self.parameterfile} does not exist.")
        if self.statefile is not None and not self.statefile.exists():
            raise FileNotFoundError(f"State file {self.statefile} does not exist.")
        self.parameters = gvec.util.read_parameters(
            self.parameterfile, format="ini"
        )  # read-only !!!

    def __del__(self):
        logger.debug(f"Deleting state {self!r}")
        if _binding.bound_state is self:
            self.unbind()
        if self._stdout is not None:
            self._stdout.close()

    # === Additional Constructors === #

    @classmethod
    def new(cls, parameters: Mapping, rundir: str | Path | None = None):
        if rundir is None:
            rundir = Path.cwd()
        parameterfile = Path(rundir) / "parameter.ini"
        gvec.util.write_parameters(parameters, parameterfile)
        return cls(parameterfile)

    # === Binding to the Fortran library === #

    from ._binding import bind, unbind

    # === Context Manager === #

    def __enter__(self):
        warnings.warn(
            "Using State as a context manager is no longer necessary.", DeprecationWarning
        )
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    # === Debug Information === #

    def __repr__(self):
        return (
            "<pygvec.State("
            + ",".join(
                [
                    str(self.parameterfile),
                    str(self.statefile),
                ]
            )
            + ")>"
        )

    @property
    def stdout(self) -> str:
        if self._stdout is None:
            return ""
        if _binding.bound_state is self:
            flush_stdout()
        self._stdout.seek(0)
        return self._stdout.read()

    @property
    def rundir(self):
        return self.parameterfile.parent

    @property
    def name(self):
        """The name of the configuration / `ProjectName` in the parameter file."""
        return self.parameters.get("ProjectName", "GVEC")

    # === Properties === #

    from ._properties import (
        nfp,
        get_integration_points,
        get_radial_gridpoints,
        get_mn_max,
    )

    # === Evaluation Methods === #

    from ._evaluate import (
        evaluate_base_tens,
        evaluate_base_tens_all,
        evaluate_base_list_tz,
        evaluate_base_list_tz_all,
        evaluate_base_list_rtz_all,
        evaluate_hmap,
        evaluate_hmap_only,
        evaluate_hmap_derivs,
        evaluate_metric_derivs,
        evaluate_jac_h_derivs,
        evaluate_profile,
        evaluate_rho2_profile,
    )

    # === Boozer Potential === #

    from ._sfl import get_boozer, get_boozer_angles, evaluate_boozer_list_tz_all

    # === Straight-Fieldline PEST angles === #

    from ._sfl import get_pest_angles

    # === High Level Interface for Evaluations === #

    def compute(
        self,
        ev: xr.Dataset,
        *quantities: str,
    ):
        """
        Compute the target equilibrium quantity and add it to the given evaluation dataset.

        This method will recursively determine prerequisites, compute them and add them to the dataset as needed.

        Parameters
        ----------
        ev : xr.Dataset
            The evaluation dataset with the target grid ``(rad, pol, tor)``, coordinates ``(rho, theta, zeta)`` and possibly some precomputed quantities.
        *quantities : str
            One or more names of the quantities to compute.
            See the :ref:`default table of available quantities <table-of-quantities>`
            or call ``table_of_quantities`` to see all options.

        See Also
        --------
        gvec.core.compute.compute: this function as a standalone function.
        evaluate: create a new grid in logical coordinates and compute target quantities.
        evaluate_sfl: create a new grid in straight-fieldline coordinates and compute target quantities.
        """
        from gvec.core.compute import compute

        return compute(ev, *quantities, state=self)

    # def evaluate(...) -> injected in gvec.core.compute
    # def evaluate_sfl(...) -> injected in gvec.core.compute
