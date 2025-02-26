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

from . import lib
from .lib import modgvec_py_post as _post
from .lib import modgvec_py_binding as _binding

from pathlib import Path
from typing import Mapping, Callable, Iterable, Literal
import re
import inspect
import functools
import tempfile
import logging
import os

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
    func: callable, argnames: Iterable[str], n_out: int, vector_out: bool = False
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
    # === Constructor & Destructor === #

    def __init__(
        self,
        parameterfile: str | Path,
        statefile: str | Path,
        redirect_stdout: bool = True,
    ):
        self.initialized: bool = False
        self.parameterfile: Path | None = None
        self.statefile: Path | None = None
        self.logger = logging.getLogger("pyGVEC.State")

        if _post.initialized:
            raise NotImplementedError("Only one instance of State is allowed.")
        if not Path(parameterfile).exists():
            raise FileNotFoundError(f"Parameter file {parameterfile} does not exist.")
        if not Path(statefile).exists():
            raise FileNotFoundError(f"State file {statefile} does not exist.")

        _binding.redirect_abort()  # redirect abort to raise a RuntimeError
        if redirect_stdout:
            self._stdout = tempfile.NamedTemporaryFile(mode="r", prefix="gvec-stdout-")
            _binding.redirect_stdout(self._stdout.name)
        self.parameterfile = Path(parameterfile)
        self.statefile = Path(statefile)

        self._original_dir = os.getcwd()
        os.chdir(self.parameterfile.parent)
        self.statefile = self.statefile.relative_to(self.parameterfile.parent)
        self.parameterfile = self.parameterfile.relative_to(self.parameterfile.parent)

        _post.init(self.parameterfile)
        _post.readstate(self.statefile)
        self.initialized = True
        self._children = []

    @_assert_init
    def finalize(self):
        """Finalize the state and free all (fortran) resources."""
        self.logger.debug(f"Finalizing state {self!r}")
        for child in self._children:
            if isinstance(child, lib.Modgvec_Sfl_Boozer.t_sfl_boozer):
                if child.initialized:
                    self.logger.debug(f"Finalizing Boozer potential {child!r}")
                    child.free()
            else:
                self.logger.error(f"Unknown child: {child!r}")

        _post.finalize()
        self.initialized = False

        os.chdir(self._original_dir)

    def __del__(self):
        self.logger.debug(f"Deleting state {self!r}")
        if hasattr(self, "_stdout"):
            self._stdout.close()
        # silently ignore non-initialized states
        if self.initialized:
            self.finalize()

    # === Context Manager === #

    @_assert_init
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.logger.debug(f"Exiting context manager for state {self!r}")
        # silently ignore non-initialized states
        if self.initialized:
            self.finalize()

    # === Debug Information === #

    def __repr__(self):
        return (
            "<pygvec.State("
            + ",".join(
                [
                    "initialized" if self.initialized else "finalized",
                    self.parameterfile.name,
                    self.statefile.name,
                ]
            )
            + ")>"
        )

    @property
    def stdout(self):
        if not hasattr(self, "_stdout"):
            return None
        self._stdout.seek(0)
        return self._stdout.read()

    # === Evaluation Methods === #

    @property
    @_assert_init
    def nfp(self):
        return _post.nfp

    @_assert_init
    def get_integration_points(self, quantity: str = "LA"):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA"]:
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'."
            )

        r_n, t_n, z_n = _post.get_integration_points_num(quantity)
        r_GP, r_w = (np.zeros(r_n, dtype=np.float64) for _ in range(2))
        t_w, z_w = _post.get_integration_points(quantity, r_GP, r_w)
        return r_GP, r_w, t_n, t_w, z_n, z_w

    @_assert_init
    def get_mn_max(self, quantity: str = "all") -> tuple[int, int]:
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA", "all"]:
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'."
            )

        if quantity == "all":
            m, n = zip(*[self.get_mn_max(q) for q in ["X1", "X2", "LA"]])
            return max(m), max(n)

        return _post.get_mn_max(quantity)

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
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'."
            )
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
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'."
            )
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
        thetazeta = np.asfortranarray(thetazeta, dtype=np.float64)
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
    def evaluate_base_list_tz_all(
        self, quantity: str, rho: np.ndarray, thetazeta: np.ndarray
    ):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA"]:
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'."
            )

        rho = np.asfortranarray(rho, dtype=np.float64)
        thetazeta = np.asfortranarray(thetazeta, dtype=np.float64)
        if rho.ndim != 1:
            raise ValueError("rho must be a 1D array.")
        if thetazeta.ndim != 2 or thetazeta.shape[0] != 2:
            raise ValueError("thetazeta must be a 2D array with shape (2, n).")
        if rho.max() > 1.0 or rho.min() < 0.0:
            raise ValueError("rho must be in the range [0, 1].")

        # Q, dQ_drho, dQ_dtheta, dQ_dzeta, dQ_drr, dQ_drt, dQ_drz, dQ_dtt, dQ_dtz, dQ_dzz
        outputs = [
            np.zeros((rho.size, thetazeta.shape[1]), dtype=np.float64, order="F")
            for _ in range(10)
        ]

        _post.evaluate_base_list_tz_all(
            rho.size, thetazeta.shape[1], rho, thetazeta, quantity, *outputs
        )
        return outputs

    @_assert_init
    def evaluate_base_tens_all(
        self, quantity: str, rho: np.ndarray, theta: np.ndarray, zeta: np.ndarray
    ):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["X1", "X2", "LA"]:
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'X1', 'X2', 'LA'."
            )

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

    # === Boozer Potential === #

    @_assert_init
    def get_boozer(
        self,
        rho: np.ndarray,
        M: int | None = None,
        N: int | None = None,
        *,
        M_nyq: int | None = None,
        N_nyq: int | None = None,
        sincos: Literal["sin", "cos", "sincos"] = "sin",
        recompute_lambda: bool = True,
    ):
        r"""
        Initialize a new Boozer potential with M poloidal and N toroidal nodes for all fluxsurfaces given by rho.

        Parameters
        ----------
        M
            Number of poloidal nodes of the Boozer potential :math:`\nu_B`. Defaults to the maximum number of nodes of the basis.
        N
            Number of toroidal nodes of the Boozer potential :math:`\nu_B`. Defaults to the maximum number of nodes of the basis.
        rho
            Array of (radius-like) flux surface labels.

        Returns
        -------
        sfl_boozer
            Straight-fieldline Boozer object (wrapped Fortran object).
        """
        # --- Defaults --- #
        M_LA, N_LA = self.get_mn_max("LA")
        _, M_nyq_LA, N_nyq_LA = _post.get_integration_points_num("LA")

        if M is None:
            M = M_LA
        if N is None:
            N = N_LA
        if M_nyq is None:
            M_nyq = max(4 * M + 1, M_nyq_LA)
        if N_nyq is None:
            N_nyq = max(4 * N + 1, N_nyq_LA)

        # --- Argument Handling --- #
        if not isinstance(M, int) or not isinstance(N, int) or M < 0 or N < 0:
            raise ValueError("M and N must be non-negative integers (or None).")
        if M < M_LA or N < N_LA:
            raise ValueError(
                f"The number of poloidal and toroidal nodes for the Boozer potential must be equal or larger to those of the original lambda: ({M}, {N}) < ({M_LA}, {N_LA})"
            )
        if (
            not isinstance(M_nyq, int)
            or not isinstance(N_nyq, int)
            or M_nyq < min(2 * M + 1, M_nyq_LA)
            or N_nyq < min(2 * N + 1, N_nyq_LA)
        ):
            raise ValueError(
                f"M_nyq and N_nyq must be integers larger than min({2 * M + 1=}, {M_nyq_LA=}) and min({2 * N + 1=}, {N_nyq_LA=}) (or None)."
            )

        rho = np.asfortranarray(rho, dtype=np.float64)
        if rho.ndim != 1 or rho.max() > 1.0 or rho.min() < 1e-4:
            raise ValueError("rho must be a 1D array in the range [1e-4, 1].")

        if sincos not in ["sin", "cos", "sincos"]:
            raise ValueError("sincos must be 'sin', 'cos', or 'sincos'.")
        sincos = {"sin": " _sin_", "cos": " _cos_", "sincos": "_sin_cos_"}[sincos]

        recompute_lambda = bool(recompute_lambda)

        # --- Create & compute Boozer potential --- #
        self.logger.debug("Initializing new Boozer potential.")
        sfl_boozer = _post.init_boozer(
            (M, N), (M_nyq, N_nyq), sincos, rho.size, rho, recompute_lambda
        )
        self._children.append(sfl_boozer)
        self.logger.debug(f"Computing Boozer potential {sfl_boozer!r}")
        _post.get_boozer(sfl_boozer)

        # ToDo: wrap sfl_boozer again to make it safer?
        return sfl_boozer

    @_assert_init
    def get_boozer_angles(
        self, sfl_boozer: lib.Modgvec_Sfl_Boozer.t_sfl_boozer, tz_list: np.ndarray
    ):
        if not isinstance(sfl_boozer, lib.Modgvec_Sfl_Boozer.t_sfl_boozer):
            raise ValueError(
                f"Boozer object {sfl_boozer!r} must be of type `t_sfl_boozer`."
            )
        if sfl_boozer not in self._children:
            raise ValueError(
                f"Boozer object {sfl_boozer!r} is not known to the state {self!r}."
            )
        if not sfl_boozer.initialized:
            raise ValueError(f"Boozer object {sfl_boozer!r} is not initialized.")

        tz_list = np.asfortranarray(tz_list, dtype=np.float64)
        if tz_list.ndim != 2 or tz_list.shape[0] != 2:
            raise ValueError("thetazeta must be a 2D array with shape (2, n).")

        tz_out = np.ndarray(
            (2, tz_list.shape[1], sfl_boozer.nrho), dtype=np.float64, order="F"
        )
        sfl_boozer.find_angles(tz_list.shape[1], tz_list, tz_out)
        return tz_out

    @_assert_init
    def evaluate_boozer_list_tz_all(
        self,
        sfl_boozer: lib.Modgvec_Sfl_Boozer.t_sfl_boozer,
        quantity: str,
        rad: np.ndarray,
        thetazeta: np.ndarray,
    ):
        if not isinstance(quantity, str):
            raise ValueError("Quantity must be a string.")
        elif quantity not in ["LA", "NU"]:
            raise ValueError(
                f"Unknown quantity: {quantity}, expected one of 'LA', 'NU'."
            )
        if not isinstance(sfl_boozer, lib.Modgvec_Sfl_Boozer.t_sfl_boozer):
            raise ValueError(
                f"Boozer object {sfl_boozer!r} must be of type `t_sfl_boozer`."
            )
        if sfl_boozer not in self._children:
            raise ValueError(
                f"Boozer object {sfl_boozer!r} is not known to the state {self!r}."
            )
        if not sfl_boozer.initialized:
            raise ValueError(f"Boozer object {sfl_boozer!r} is not initialized.")

        rad = np.asfortranarray(rad, dtype=np.int64)
        thetazeta = np.asfortranarray(thetazeta, dtype=np.float64)
        if rad.ndim != 1:
            raise ValueError("rad must be a 1D array.")
        if thetazeta.ndim != 2 or thetazeta.shape[0] != 2:
            raise ValueError("thetazeta must be a 2D array with shape (2, n).")
        if rad.min() < 0:
            raise ValueError("rad must be a positive integer.")

        # Q, dQ_dtheta, dQ_dzeta, dQ_dtt, dQ_dtz, dQ_dzz
        outputs = [
            np.zeros((rad.size, thetazeta.shape[1]), dtype=np.float64, order="F")
            for _ in range(6)
        ]

        _post.evaluate_boozer_list_tz_all(
            sfl_boozer, rad.size, thetazeta.shape[1], rad, thetazeta, quantity, *outputs
        )
        return outputs

    # === Integration with computable quantities === #

    def compute(self, ds: xr.Dataset, *quantities):
        from .comp import compute

        return compute(ds, *quantities, state=self)
