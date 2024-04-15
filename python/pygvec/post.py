# ============================================================================================================================== #
# Copyright (C) 2024 Robert Babin <robert.babin@ipp.mpg.de>
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
"""pygvec postprocessing"""

from ._post import modpygvec_post as _post

from pathlib import Path
import numpy as np
from typing import Sequence
from collections import Counter

_status = "pre"


class State:
    def __init__(self, parameterfile, statefile):
        global _status
        if _status != "pre":
            raise NotImplementedError("Only one instance of State is allowed.")
        _status = "init"
        if not Path(parameterfile).exists():
            raise FileNotFoundError(f"Parameter file {parameterfile} does not exist.")
        if not Path(statefile).exists():
            raise FileNotFoundError(f"State file {statefile} does not exist.")

        _post.init(parameterfile)
        _post.readstate(statefile)

    def __enter__(self):
        global _status
        if _status != "init":
            raise NotImplementedError("State is not initialized.")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        global _status
        if _status == "init":
            _post.finalize()
            _status = "final"

    def __del__(self):
        global _status
        if _status == "init":
            _post.Finalize()
            _status = "final"

    def evaluate_base(
        self,
        rho: np.ndarray | Sequence,
        theta: np.ndarray | Sequence,
        zeta: np.ndarray | Sequence,
        quantity: str,
    ):
        global _status
        if _status != "init":
            raise NotImplementedError("State is not initialized.")

        if quantity in ["X1", "X2", "LA"]:
            selection = (quantity, "", "")
        elif quantity.startswith("D_"):
            deriv, quantity2 = quantity[2:].split(" ")
            if quantity2 not in ["X1", "X2", "LA"]:
                raise ValueError(f"Unknown quantity: {quantity2} in {quantity}")
            counts = Counter(deriv)
            if set(counts.keys()) - set("rtz"):
                raise ValueError(f"Unknown derivative: {deriv} in {quantity}")
            if (
                counts["r"] > 2
                or counts["t"] > 2
                or counts["z"] > 2
                or counts["t"] + counts["z"] > 2
            ):
                raise ValueError(
                    f"Derivative {deriv} is too high. (max. 'rr' and 'tt', 'tz', 'zz' respectively)"
                )
            selection = (
                quantity2,
                "s" * counts["r"],
                "t" * counts["t"] + "z" * counts["z"],
            )
        else:
            raise ValueError(f"Unknown quantity: {quantity}")

        
        rho = np.asarray(rho, dtype=np.float64) # Fortran order does not matter (1D)
        theta = np.asarray(theta, dtype=np.float64)
        zeta = np.asarray(zeta, dtype=np.float64)
        if rho.ndim != 1 or theta.ndim != 1 or zeta.ndim != 1:
            raise ValueError("rho, theta, and zeta must be 1D arrays.")
        if rho.max() > 1.0 or rho.min() < 0.0:
            raise ValueError("rho must be in the range [0, 1].")

        result = np.zeros(
            (zeta.size, theta.size, rho.size), dtype=np.float64, order="F"
        )
        _post.evaluate_base(rho, theta, zeta, *selection, result)
        return result.T  # C order (rho, theta, zeta)
