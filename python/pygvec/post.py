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
from collections import Counter
from typing import Sequence, Mapping
import re
import inspect
import functools

import numpy as np
import xarray as xr


def _assert_init(func):
    @functools.wraps(func)
    def wrapped(self, *args, **kwargs):
        if not self.initialized:
            raise RuntimeError("State is not initialized.")
        if not _post.initialized:
            raise RuntimeError("State is initialized, but GVEC libaray is not!")
        return func(self, *args, **kwargs)

    return wrapped


def _evaluate_1D_factory(
    func: callable, argnames: Sequence[str], n_out: int, vector_out: bool = False
):
    params = [inspect.Parameter("self", inspect.Parameter.POSITIONAL_OR_KEYWORD)] + [
        inspect.Parameter(
            name, inspect.Parameter.POSITIONAL_OR_KEYWORD, annotation=np.ndarray
        )
        for name in argnames
    ]
    returns = tuple[tuple(np.ndarray for _ in range(n_out))]
    sig = inspect.Signature(params, return_annotation=returns)

    @_assert_init
    def wrapper(self, *args, **kwargs):
        bound_args = sig.bind(self, *args, **kwargs)
        inputs = [
            np.asfortranarray(value, dtype=np.float64)
            for key, value in bound_args.arguments.items()
            if key != "self"
        ]
        n = inputs[0].size
        for value in inputs:
            if value.shape != (n,):
                raise ValueError("All arguments must be 1D arrays of the same size.")

        if vector_out:
            outputs = [
                np.zeros((3, n), dtype=np.float64, order="F") for _ in range(n_out)
            ]
        else:
            outputs = [np.zeros(n, dtype=np.float64) for _ in range(n_out)]
        func(n, *inputs, *outputs)
        return outputs

    wrapper.__signature__ = sig
    wrapper.__name__ = func.__name__
    return wrapper


class State:
    def __init__(self, parameterfile, statefile):
        self.initialized = False
        self.parameterfile = None
        self.statefile = None

        if _post.initialized:
            raise NotImplementedError("Only one instance of State is allowed.")
        if not Path(parameterfile).exists():
            raise FileNotFoundError(f"Parameter file {parameterfile} does not exist.")
        if not Path(statefile).exists():
            raise FileNotFoundError(f"State file {statefile} does not exist.")

        self.parameterfile = Path(parameterfile)
        _post.init(parameterfile)
        self.statefile = Path(statefile)
        _post.readstate(statefile)
        self.initialized = True

    @_assert_init
    def finalize(self):
        _post.finalize()
        self.initialized = False

    @_assert_init
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # silently ignore non-initialized states
        if self.initialized:
            self.finalize()

    def __del__(self):
        # silently ignore non-initialized states
        if self.initialized:
            self.finalize()

    def __repr__(self):
        return f"<pygvec.State({bool(self.initialized)},{self.parameterfile},{self.statefile})>"

    @_assert_init
    def evaluate_base_tens(
        self,
        quantity: str,
        derivs: str | None,
        rho: np.ndarray,
        theta: np.ndarray,
        zeta: np.ndarray,
    ):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA"]:
            raise ValueError(f"Unknown quantity: {quantity}")
        if derivs is not None:
            if not isinstance(derivs, str):
                raise ValueError("Derivatives must be a string.")
            if m := re.match(r"(r{0,2})(t{0,2}|z{0,2}|tz)$", derivs):
                sel_derivs = m.groups()
            else:
                raise ValueError(f"Unknown derivative: {derivs}")
        else:
            sel_derivs = ("", "")

        rho = np.asfortranarray(rho, dtype=np.float64)
        theta = np.asfortranarray(theta, dtype=np.float64)
        zeta = np.asfortranarray(zeta, dtype=np.float64)
        if rho.ndim != 1 or theta.ndim != 1 or zeta.ndim != 1:
            raise ValueError("rho, theta, and zeta must be 1D arrays.")
        if rho.max() > 1.0 or rho.min() < 0.0:
            raise ValueError("rho must be in the range [0, 1].")

        result = np.zeros(
            (rho.size, theta.size, zeta.size), dtype=np.float64, order="F"
        )
        _post.evaluate_base_tens(rho, theta, zeta, quantity, *sel_derivs, result)
        return result

    @_assert_init
    def evaluate_base_list_tz(
        self, quantity: str, derivs: str | None, rho: np.ndarray, thetazeta: np.ndarray
    ):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA"]:
            raise ValueError(f"Unknown quantity: {quantity}")
        if derivs is not None:
            if not isinstance(derivs, str):
                raise ValueError("Derivatives must be a string.")
            if m := re.match(r"(r{0,2})(t{0,2}|z{0,2}|tz)$", derivs):
                sel_derivs = m.groups()
            else:
                raise ValueError(f"Unknown derivative: {derivs}")
        else:
            sel_derivs = ("", "")

        rho = np.asfortranarray(rho, dtype=np.float64)
        thetazeta = np.asfortranarray(thetazeta, dtype=np.float64, order="F")
        if rho.ndim != 1:
            raise ValueError("rho must be a 1D array.")
        if thetazeta.ndim != 2 or thetazeta.shape[0] != 2:
            raise ValueError("thetazeta must be a 2D array with shape (2, n).")
        if rho.max() > 1.0 or rho.min() < 0.0:
            raise ValueError("rho must be in the range [0, 1].")

        result = np.zeros((rho.size, thetazeta.shape[1]), dtype=np.float64, order="F")
        _post.evaluate_base_list_tz(
            rho.size, thetazeta.shape[1], rho, thetazeta, quantity, *sel_derivs, result
        )
        return result

    @_assert_init
    def evaluate_base_tens_all(
        self, quantity: str, rho: np.ndarray, theta: np.ndarray, zeta: np.ndarray
    ):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA"]:
            raise ValueError(f"Unknown quantity: {quantity}")

        rho = np.asfortranarray(rho, dtype=np.float64)
        theta = np.asfortranarray(theta, dtype=np.float64)
        zeta = np.asfortranarray(zeta, dtype=np.float64)
        if rho.ndim != 1 or theta.ndim != 1 or zeta.ndim != 1:
            raise ValueError("rho, theta, and zeta must be 1D arrays.")
        if rho.max() > 1.0 or rho.min() < 0.0:
            raise ValueError("rho must be in the range [0, 1].")

        # Q, dQ_drho, dQ_dtheta, dQ_dzeta, dQ_drr, dQ_drt, dQ_drz, dQ_dtt, dQ_dtz, dQ_dzz
        outputs = [
            np.zeros((rho.size, theta.size, zeta.size), dtype=np.float64, order="F")
            for _ in range(10)
        ]

        _post.evaluate_base_tens_all(
            rho.size, theta.size, zeta.size, rho, theta, zeta, quantity, *outputs
        )
        return outputs

    evaluate_hmap = _evaluate_1D_factory(
        _post.evaluate_hmap,
        [
            "X1",
            "X2",
            "zeta",
        ]
        + [f"d{Q}_d{i}" for i in "rtz" for Q in ["X1", "X2"]],
        4,
        True,
    )  # -> pos, e_rho, e_theta, e_zeta

    evaluate_hmap_only = _evaluate_1D_factory(
        _post.evaluate_hmap_only, ["X1", "X2", "zeta"], 4, True
    )  # -> pos, e_X1, e_X2, e_zeta3

    evaluate_metric = _evaluate_1D_factory(
        _post.evaluate_metric,
        [
            "X1",
            "X2",
            "zeta",
        ]
        + [
            f"d{Q}_d{i}"
            for i in "r t z rr rt rz tt tz zz".split()
            for Q in ["X1", "X2"]
        ],
        24,
    )  # -> g_rr, g_rt ... g_zz, dg_rr_dr, dg_rt_dr ... dg_zz_dz

    evaluate_jacobian = _evaluate_1D_factory(
        _post.evaluate_jacobian,
        [
            "X1",
            "X2",
            "zeta",
            "dX1_dr",
            "dX2_dr",
            "dX1_dt",
            "dX2_dt",
            "dX1_dz",
            "dX2_dz",
        ],
        4,
    )  # -> Jac_h, dJac_h_dr, dJac_h_dt, dJac_h_dz

    @_assert_init
    def evaluate_profile(self, quantity: str, rho: np.ndarray):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in [
            "iota",
            "iota_prime",
            "p",
            "p_prime",
            "chi",
            "chi_prime",
            "Phi",
            "Phi_prime",
            "Phi_2prime",
            "PhiNorm",
            "PhiNorm_prime",
        ]:
            raise ValueError(f"Unknown quantity: {quantity}")

        rho = np.asfortranarray(rho, dtype=np.float64)
        if rho.ndim != 1:
            raise ValueError("rho must be a 1D array.")
        if rho.max() > 1.0 or rho.min() < 0.0:
            raise ValueError("rho must be in the range [0, 1].")

        result = np.zeros(rho.size, dtype=np.float64, order="F")
        _post.evaluate_profile(rho.size, rho, quantity, result)
        return result


class Evaluations(xr.Dataset):
    __slots__ = ["state"]
    _quantities = {}

    def __init__(self, state: State, **kwargs):
        self.state = state
        super().__init__(**kwargs)

    @classmethod
    def register_quantity(
        cls,
        name: None | str | Sequence[str] = None,
        coords: Sequence[str] = (...,),
        requirements: Sequence[str] = (),
    ):
        """Decorator to register equilibrium quantities.

        The decorated functions should have the signature `func(xarray.Dataset, State) -> xarray.Dataset`.
        """

        def _register(func):
            nonlocal name
            if name is None:
                name = func.__name__
            func.name = name
            func.coords = coords
            func.requirements = requirements
            if isinstance(name, str):
                name = [name]
            for n in name:
                if n in cls._quantities:
                    raise KeyError(
                        f"A quantity `{n}` is already registered with {cls}."
                    )
                cls._quantities[n] = func
            return func

        return _register

    def compute(self, quantity: str):
        """Compute the target equilibrium quantity.

        This method will compute required parameters recursively.
        """
        if quantity in self:
            return self
        if quantity not in self._quantities:
            raise KeyError(f"The quantity `{quantity}` is not registered.")
        func = self._quantities[quantity]
        if "vector" in func.coords and "vector" not in self.coords:
            self.coords["vector"] = ["x", "y", "z"]
        for coord in func.coords:
            if coord != ... and coord not in self.coords:
                raise ValueError(
                    f"The quantity `{quantity}` requires a grid of `({', '.join(map(str, func.coords))})` coordinates."
                )
        for req in func.requirements:
            if req not in self:
                if req == "mu0":
                    self["mu0"] = 4 * np.pi * 1e-7
                    self.mu0.attrs["long_name"] = "magnetic constant"
                    self.mu0.attrs["symbol"] = r"\mu_0"
                self.compute(req)

        if "state" in inspect.signature(func).parameters:
            func(self, self.state)
        else:
            func(self)
