import numpy as np
from gvec.core.compute import fluxsurface_integral
from gvec.core.run_state import RunState
from gvec.core.state import with_binding
from typing import Callable

from gvec.lib import modgvec_py_state as _state

import logging
from pathlib import Path
from scipy.constants import mu_0


logger = logging.getLogger("gvec.free-boundary")


class FreeBoundaryState(RunState):
    def __init__(
        self,
        parameterfile: str | Path,
        statefile: str | Path | None,
        I_tor: float | Callable | None,
        reset_LA: bool,
        n_iota_rhos: int,
        stdout_path: str | None,
        dt: float,
    ):
        """Free-boundary state object.

        Parameters
        ----------
        parameterfile : str | Path
            Parameter.ini file used for initialization.
        statefile : str | Path | None, optional
            Statefile to restart from, by default None
        """
        super().__init__(parameterfile, statefile, I_tor, reset_LA, n_iota_rhos, stdout_path)

        # initialize diagnostics dataset for quantities at the boundary for each free-boundary update

        r_GP, r_w, t_n, t_w, z_n, z_w = self.get_integration_points()

        self.t_n = t_n
        self.z_n = z_n

        self.dt = dt

        self.m_max, self.n_max = self.get_mn_max()
        self.clockwise_zeta = self._determine_clockwise_zeta()

    @with_binding
    def set_Bext_squared(self, B2):
        logger.debug("Updating exterior field")
        _state.set_bext_squared(B2)

    @with_binding
    def _get_boundary_forces(self):
        n_X1 = self._get_X1_n_modes()
        n_X2 = self._get_X2_n_modes()
        forces_out = np.asfortranarray(np.zeros((3, n_X1 + n_X2)))
        _state.get_boundary_forces(n_X1, n_X2, forces_out)
        m = forces_out[0, :]
        n = forces_out[1, :]
        F = forces_out[2, :]
        return (m, n, F)

    @with_binding
    def get_Fnorm_boundary(self):
        logger.debug("Extracting Fnorm at the boundary")
        return _state.get_fnorm_boundary()

    @with_binding
    def get_gnorm(self):
        logger.debug("Extracting gnorm")
        gnorm = _state.get_gnorm()
        return gnorm

    @with_binding
    def get_xnorm(self):
        logger.debug("Extracting xnorm")
        xnorm = _state.get_xnorm()
        return xnorm

    @with_binding
    def update_boundary(self, dt: float, reset_minimizer: bool, blend_depth: float = 0.6):
        logger.debug("Stepping in free boundary direction")
        dt = _state.update_boundary(dt, reset_minimizer, blend_depth=blend_depth)
        return dt

    def _determine_clockwise_zeta(self):
        """Rotates a point of the axis with the Rodrigues equation to see if zeta runs clockwise or counter-clockwise.

        Returns
        -------
        bool
            Whether zeta runs clockwise ort not.

        Raises
        ------
        ValueError
            If there is an implementation error that causes this approach to fail.
        """

        zeta = 2 * np.pi / self.nfp
        ev = self.evaluate("pos", rho=0.0, theta=0.0, zeta=[0, 2 * np.pi / self.nfp])
        v0 = ev.pos.isel(pol=0, tor=0, rad=0).data
        v1 = ev.pos.isel(pol=0, tor=1, rad=0).data
        k = np.array([0, 0, 1])
        vrot = (
            v0 * np.cos(zeta)
            + np.sin(zeta) * np.cross(k, v0)
            + (1 - np.cos(zeta)) * np.dot(k, v0) * k
        )
        zeta_inv = 2 * np.pi - zeta
        vrot_inv = (
            v0 * np.cos(zeta_inv)
            + np.sin(zeta_inv) * np.cross(k, v0)
            + (1 - np.cos(zeta_inv)) * np.dot(k, v0) * k
        )

        if all(abs(vrot - v1) < 1e-12):
            clockwise_zeta = False
        elif all(abs(vrot_inv - v1) < 1e-12):
            clockwise_zeta = True
        else:
            raise ValueError("Issue when determining zeta-rotation.")
        return clockwise_zeta

    def _free_boundary_target(self, B_ext_2: np.ndarray):
        """Evaluate the free-boundary target function based on the force-balance.

        Parameters
        ----------
        B_ext_2 : np.ndarray
            Squared external magnetic field

        Returns
        -------
        float
            free-boundary target
        """
        ev_gvec = (
            self.evaluate("B", "p", "normal", "Phi_edge", rho=1.0)
            .sel(rho=1.0)
            .transpose("xyz", "tor", "pol")
        )
        B_gvec_2 = (ev_gvec.B**2).sum(dim="xyz")
        force_balance = (B_gvec_2 - B_ext_2) / 2 + mu_0 * ev_gvec.p

        B_avg = (B_gvec_2 + B_ext_2) / 2
        # flux surface averaged normalized force balance
        # <|[[B²/2+mu_0*P]] / B² |>
        free_boundary_target = fluxsurface_integral(
            abs(force_balance / B_avg) * ev_gvec.Jac
        ) / fluxsurface_integral(ev_gvec.Jac)
        return free_boundary_target, force_balance, ev_gvec

    @with_binding
    def _free_boundary_update(
        self,
        maxIter: int,
        B2: np.ndarray,
        reset_minimizer: bool = False,
        blend_depth: float = 0.6,
        **minimize_kwargs,
    ):
        """Single free boundary push.

        Updates the FORTRAN B_ext^2, resets the fixed-boundary minimizer,
        pushes the boundary and minimizes the fixed-boundary problem

        Parameters
        ----------
        maxIter : int
            Maximum number of iterations for the fixed-boundary minimization.
        B2 : np.ndarray
            External magnetic field, flattened in GVEC order (zeta->theta)
        tol : float
            Tolerance on fixed boundary forces.
        copy_stdout : bool, optional
            Whether to save stdout to stoudt.txt, by default True
        reset_minimizer: bool, optional
            Whether to reset the free boundary minimizers to their initial states.
        """
        # set Bext
        self.set_Bext_squared(B2)
        self.reset_descent()
        self.dt = self.update_boundary(self.dt, reset_minimizer, blend_depth)

        self.minimize(logfinaliter=False, maxIter=maxIter, **minimize_kwargs)
