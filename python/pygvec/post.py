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
from collections import Counter
import xarray as xr

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
            _post.finalize()
            _status = "final"

    def evaluate_base(
        self,
        quantity: str,
        *args,
        **kwargs,
    ):
        global _status
        if _status != "init":
            raise NotImplementedError("State is not initialized.")

        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity in ["all", "X1", "X2", "LA"]:
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

        if len(args) == 1 and len(kwargs) == 0:
            if not isinstance(args[0], xr.Dataset):
                raise ValueError("'ds' must be an xarray.DataArray.")
            ds = args[0]
            rho = ds.coords["rho"].values
            theta = ds.coords["theta"].values
            zeta = ds.coords["zeta"].values
            algorithm = "tensorproduct"
        elif len(args) == 0 and {"ds"} == set(kwargs.keys()):
            if not isinstance(kwargs["ds"], xr.DataArray):
                raise ValueError("'ds' must be an xarray.DataArray.")
            ds = kwargs["ds"]
            rho = ds.coords["rho"].values
            theta = ds.coords["theta"].values
            zeta = ds.coords["zeta"].values
            algorithm = "tensorproduct"
        elif len(args) == 3 and len(kwargs) == 0:
            rho, theta, zeta = args
            algorithm = "tensorproduct"
        elif len(args) == 0 and {"rho", "theta", "zeta"} == set(kwargs.keys()):
            rho = kwargs["rho"]
            theta = kwargs["theta"]
            zeta = kwargs["zeta"]
            algorithm = "tensorproduct"
        elif len(args) == 0 and {"r", "t", "z"} == set(kwargs.keys()):
            rho = kwargs["r"]
            theta = kwargs["t"]
            zeta = kwargs["z"]
            algorithm = "tensorproduct"
        elif len(args) == 2 and len(kwargs) == 0:
            rho, thetazeta = args
            algorithm = "list-tz"
        elif len(args) == 0 and {"rho", "thetazeta"} == set(kwargs.keys()):
            rho = kwargs["rho"]
            thetazeta = kwargs["thetazeta"]
            algorithm = "list-tz"
        elif len(args) == 0 and {"r", "tz"} == set(kwargs.keys()):
            rho = kwargs["r"]
            thetazeta = kwargs["tz"]
            algorithm = "list-tz"
        else:
            raise ValueError("Invalid arguments.")

        if algorithm == "tensorproduct":
            rho = np.asfortranarray(rho, dtype=np.float64)
            theta = np.asfortranarray(theta, dtype=np.float64)
            zeta = np.asfortranarray(zeta, dtype=np.float64)
            if rho.ndim != 1 or theta.ndim != 1 or zeta.ndim != 1:
                raise ValueError("rho, theta, and zeta must be 1D arrays.")
            if rho.max() > 1.0 or rho.min() < 0.0:
                raise ValueError("rho must be in the range [0, 1].")

            if quantity == "all":
                # X1, X2, dX1_ds, dX2_ds, dX1_dthet, dX2_dthet, dX1_dzeta, dX2_dzeta
                results = [
                    np.zeros(
                        (rho.size, theta.size, zeta.size), dtype=np.float64, order="F"
                    )
                    for _ in range(8)
                ]
                _post.evaluate_base_tens_all(
                    rho.size, theta.size, zeta.size, rho, theta, zeta, *results
                )
                if "ds" in locals():
                    return ds.assign(
                        {
                            "X1": (["rho", "theta", "zeta"], results[0]),
                            "X2": (["rho", "theta", "zeta"], results[1]),
                            "dX1_drho": (["rho", "theta", "zeta"], results[2]),
                            "dX2_drho": (["rho", "theta", "zeta"], results[3]),
                            "dX1_dtheta": (["rho", "theta", "zeta"], results[4]),
                            "dX2_dtheta": (["rho", "theta", "zeta"], results[5]),
                            "dX1_dzeta": (["rho", "theta", "zeta"], results[6]),
                            "dX2_dzeta": (["rho", "theta", "zeta"], results[7]),
                        }
                    )
                return results
            else:
                result = np.zeros(
                    (rho.size, theta.size, zeta.size), dtype=np.float64, order="F"
                )
                _post.evaluate_base_tens(rho, theta, zeta, *selection, result)
                if "ds" in locals():
                    return ds.assign({quantity: (["rho", "theta", "zeta"], result)})
                return result
        elif algorithm == "list-tz":
            rho = np.asfortranarray(rho, dtype=np.float64)
            thetazeta = np.asfortranarray(thetazeta, dtype=np.float64, order="F")
            if rho.ndim != 1:
                raise ValueError("rho must be a 1D array.")
            if thetazeta.ndim != 2 or thetazeta.shape[0] != 2:
                raise ValueError("thetazeta must be a 2D array with shape (2, n).")
            if rho.max() > 1.0 or rho.min() < 0.0:
                raise ValueError("rho must be in the range [0, 1].")
            if quantity == "all":
                raise ValueError("Quantity 'all' is not supported for list-tz inputs.")

            result = np.zeros(
                (rho.size, thetazeta.shape[1]), dtype=np.float64, order="F"
            )
            _post.evaluate_base_list_tz(
                rho.size, thetazeta.shape[1], rho, thetazeta, *selection, result
            )
            return result
        raise RuntimeError("Unknown `algorithm`.")

    def evaluate_hmap(self, *args):
        global _status
        if _status != "init":
            raise NotImplementedError("State is not initialized.")

        # X1, X2, zeta, dX1_drho, dX2_drho, dX1_dtheta, dX2_dtheta, dX1_dzeta, dX2_dzeta
        if len(args) == 9:
            inputs = list(args)
        elif len(args) == 1:
            ds = args[0]
            if not isinstance(ds, xr.Dataset):
                raise ValueError("'ds' must be an xarray.Dataset.")
            # use .data to avoid copying, conversion to numpy is done later
            inputs = [ds[key].data for key in ["X1", "X2"]]
            inputs += [ds.zeta.broadcast_like(ds.X1).data]
            inputs += [
                ds[key].data
                for key in [
                    "dX1_drho",
                    "dX2_drho",
                    "dX1_dtheta",
                    "dX2_dtheta",
                    "dX1_dzeta",
                    "dX2_dzeta",
                ]
            ]
        else:
            raise ValueError("Invalid arguments.")
        n = inputs[0].size
        for i in range(len(inputs)):
            # assure that the array is contiguous (Fortran order isn't necessary for 1D)
            inputs[i] = np.asfortranarray(inputs[i], dtype=np.float64).flatten()
            if inputs[i].size != n:
                raise ValueError("All inputs must have the same size.")
        # coords, e_rho, e_theta, e_zeta
        outputs = [np.zeros((3, n), dtype=np.float64, order="F") for _ in range(4)]

        _post.evaluate_hmap(n, *inputs, *outputs)
        if "ds" in locals():
            ds = ds.assign_coords(vector=np.array(["x", "y", "z"]))
            return ds.assign(
                {
                    key: (
                        ["vector", "rho", "theta", "zeta"],
                        value.reshape(3, ds.rho.size, ds.theta.size, ds.zeta.size),
                    )
                    for key, value in zip(
                        ["coords", "e_rho", "e_theta", "e_zeta"], outputs
                    )
                }
            )
        return outputs
