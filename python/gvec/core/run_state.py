import numpy as np
from gvec.core.state import State
from gvec.util import read_parameters, write_parameters, CaseInsensitiveDict, chdir
from gvec.core.state import with_binding
from typing import Callable, Literal

from gvec.lib import modgvec_py_state as _state

import time
import logging
from pathlib import Path

logger = logging.getLogger("gvec.run")


class RunState(State):
    def __init__(
        self,
        parameterfile: str | Path,
        statefile: str | Path | None,
        I_tor: float | Callable | None,
        reset_LA: bool,
        n_iota_rhos: int,
        stdout_path: str | None,
    ):
        super().__init__(parameterfile, statefile, init_for_minimization=True)
        self._stdout = stdout_path

        # these parameters are modified during current constraint!
        self._parameter_dict = read_parameters(parameterfile)

        self.reset_LA = reset_LA

        if I_tor is not None:
            self.current_constraint = True
            self.rho2_iota = np.linspace(0, 1, n_iota_rhos) ** 2
            if callable(I_tor):
                self.I_tor = I_tor(self.rho2_iota)
            else:
                self.I_tor = I_tor * np.ones_like(self.rho2_iota)
            self.I_tor_boundary = self.I_tor[-1]
        else:
            self.current_constraint = False
            self.rho2_iota = None
            self.I_tor = None

        # initialize diagnostics dataset for quantities at the boundary for each free-boundary update

        r_GP, r_w, t_n, t_w, z_n, z_w = self.get_integration_points()

        self.t_n = t_n
        self.z_n = z_n

        self.iters_used = 0
        self.rms_iota = None
        self.delta_iota = None

    @with_binding
    def _minimize(self, dt=None, logfinaliter=True, maxIter=None, minimize_tol=None):
        """Minimize the MHD functional for the current state."""
        logger.debug("Minimizing MHD functional.")
        with chdir(self.rundir):
            _state.minimize(dt, logfinaliter, maxIter, minimize_tol)

    @with_binding
    def save_statefile(self):
        logger.debug("Writing state to file")
        _state.save_statefile()

    @with_binding
    def update_iota_for_current_constraint(self, rho2, iota, reset_LA):
        logger.debug("Updating iota")
        _state.update_iota_for_current_constraint(rho2, iota, reset_LA)

    @with_binding
    def get_outputlevel(self):
        logger.debug("Extracting outputlevel")
        outid = _state.get_outputlevel()
        return outid

    @with_binding
    def get_iterations(self):
        logger.debug("Extracting iterations")
        iters = _state.get_iterations()
        return iters

    @with_binding
    def get_minimize_tol(self):
        logger.debug("Extracting minimize tol")
        tol = _state.get_minimize_tol()
        return tol

    @with_binding
    def reset_descent(self):
        logger.debug("Resetting accelerated descent")
        _state.reset_descent()

    @with_binding
    def save(self):
        self.save_statefile()
        iters = self.get_iterations()
        outid = self.get_outputlevel()
        write_parameters(
            self._parameter_dict,
            self._parameter_dict["ProjectName"] + f"_parameters_{outid:04}_{iters:08}.ini",
        )

    @with_binding
    def _get_new_iota(self):
        """Calculate a new iota profile. Extracted from gvec-run, to be used for current constraint

        Returns
        -------
        (xr.DataArray, xr.DataArray, xr.DataArray)
            new iota values, rho^2 interpolation points, root mean square iota error
        """
        rho_eval = np.concatenate(
            [[np.sqrt(1e-8), np.sqrt(2e-8), np.sqrt(3e-8)], np.sqrt(self.rho2_iota[1:])]
        )
        quantites = ["iota", "iota_curr_0", "iota_0", "I_tor"]
        ev = self.evaluate(*quantites, rho=rho_eval)
        ev = ev[quantites]
        r1 = ev.isel(rad=0)
        r2 = ev.isel(rad=1)
        r3 = ev.isel(rad=2)
        ev = ev.isel(rad=slice(2, None))

        # workaround to modify xarray coordinate & index (with pandas >=3.0)
        rho = ev.rho.data.copy()
        rho[0] = 0.0
        ev = ev.assign_coords(rho=("rad", rho)).set_xindex("rho")

        # extrapolation to axis
        for var in ev.data_vars:
            ev[var].data[0] = 3 * (r1[var].data - r2[var].data) + r3[var].data

        iota_values = ev.iota_0 + self.I_tor * ev.iota_curr_0
        delta_iota = ev.iota - iota_values
        rms_iota = np.sqrt((delta_iota**2).mean(dim="rad")).item()
        return (iota_values, ev.rho**2, rms_iota, delta_iota)

    @with_binding
    def _minimize_single_current_constraint(
        self,
        maxIter: int,
        maxIter_per_update: int,
        iota_tol: float,
        minimize_tol: float | None,
        maxUpdates: int,
        reset_LA: bool,
    ):
        """Single current optimization stage.

        Parameters
        ----------
        maxIter : int
            Maximum number of fixed-boundary iterations for this stage.
        maxIter_per_update : int, optional
            Maximum number of fixed-boundary iterations per iota update, by default int(1e12).
        iota_tol : float, optional
            Tolerance on the RMS of the iota update, by default 1e-10.
        minimize_tol : float | None, optional
            Tolerance for the norm of the forces, by default None
        maxUpdates : int, optional
            Maximum number of iota updates

        Returns
        -------
        float
            RMS of the iota update.
        """
        rms_iota = 1e6
        nth_update = -1
        iters_used_start = self.get_iterations()
        self.reset_descent()
        iters_used = 0
        delta_iota = 1e6
        delta_iota_over_run = []
        iters_per_run = []

        while (rms_iota > iota_tol) and (iters_used < maxIter):
            if nth_update + 1 > maxUpdates:
                logger.warning("Maximum number of iota updates reached reached!")
                break
            self._minimize(
                logfinaliter=False,
                maxIter=min(maxIter - iters_used, maxIter_per_update),
                minimize_tol=minimize_tol,
            )
            new_iota, new_rho2, rms_iota, delta_iota = self._get_new_iota()

            # only reset LA at the start
            if iters_used > 0:
                reset_LA = False

            self.update_iota_for_current_constraint(new_rho2, new_iota, reset_LA)
            self._parameter_dict["iota"] = CaseInsensitiveDict(
                type="interpolation",
                rho2=new_rho2.data,
                vals=new_iota.data,
            )
            self.iota = new_iota

            iters = self.get_iterations()
            iters_used = iters - iters_used_start

            nth_update += 1
            delta_iota_over_run.append(rms_iota)
            iters_per_run.append(iters)

        if rms_iota > iota_tol and maxIter > 0:
            logger.warning(
                "Targeted iota has not been reached!\n"
                + f"iota_tol.: {iota_tol:.2e}, achieved rms Δiota.: {rms_iota:.2e}"
            )
        logger.info(
            f"Terminated current minimization after {nth_update} updates. Δiota = {rms_iota:.3e}, target Δiota = {iota_tol:.3e}. Used {iters_used} / {maxIter} iterations."
        )
        return (
            rms_iota,
            delta_iota,
            np.array(delta_iota_over_run),
            np.array(iters_per_run),
            nth_update,
        )

    @with_binding
    def _minimize_current_constraint(
        self,
        iota_target: float,
        minimize_tol: float,
        maxIter: int,
        maxUpdates: int,
        reset_LA: bool,
        target: Literal["iota", "iota_and_force", "auto"],
        **kwargs,
    ):
        """Fixed-boundary run with prescribed current profile.

        Parameters
        ----------
        iota_target : float
            Tolerance on the RMS of the iota update.
        minimize_tol : float
            Tolerance for the norm of the forces.
        maxIter : int
            Maximum number of iterations for the fixed boundary run.
        maxUpdates : int, optional
            Maximum number of iota-updates per stage
        """
        if target not in ["iota", "iota_and_force", "auto"]:
            raise ValueError(
                f"""Unknown picard_current["target"]: {target}! target has to be "iota","iota_and_force" or "auto"."""
            )
        t0 = time.time()
        log_minimize_target = np.log10(minimize_tol)
        log_iota_target = np.log10(iota_target)
        n_stages = max(int(max(-2 - log_minimize_target, log_minimize_target)), 1)
        minimize_tols = np.logspace(max(-4, log_minimize_target), log_minimize_target, n_stages)
        if n_stages > 1:
            iota_tols = np.logspace(max(-3, log_iota_target), log_iota_target, n_stages)
        else:
            iota_tols = [1e-10]

        iota_tol = iota_tols[0]

        iters_used_start = self.get_iterations()

        self.nth_update = 0
        # target iota only
        if target == "iota" or target == "auto":
            (
                self.rms_iota,
                self.delta_iota,
                self.delta_iota_over_run,
                self.iters_per_run,
                nth_update,
            ) = self._minimize_single_current_constraint(
                maxIter=maxIter,
                maxIter_per_update=min(10, maxIter),
                iota_tol=iota_tol,
                minimize_tol=minimize_tols[0],
                maxUpdates=maxUpdates,
                reset_LA=reset_LA,
            )
            self.nth_update += nth_update
        else:
            self.delta_iota_over_run = np.array([])
            self.iters_per_run = np.array([])
        iters_used = self.get_iterations() - iters_used_start

        # target iota and forces
        if target == "iota_and_force" or target == "auto":
            if iters_used < maxIter:
                for minimize_tol, iota_tol in zip(minimize_tols, iota_tols):
                    logger.debug(f"minimize_tol:{minimize_tol}")
                    logger.debug(f"iota_tol:{minimize_tol}")
                    logger.debug(f"iters_used_start:{iters_used_start}")
                    logger.debug(f"iters_used:{iters_used}")
                    logger.debug(f"maxIter:{maxIter - iters_used}")
                    (
                        self.rms_iota,
                        self.delta_iota,
                        delta_iota_over_run,
                        iters_per_run,
                        nth_update,
                    ) = self._minimize_single_current_constraint(
                        maxIter=maxIter - iters_used,
                        maxIter_per_update=1e12,
                        iota_tol=iota_tol,
                        minimize_tol=minimize_tol,
                        maxUpdates=maxUpdates,
                        reset_LA=reset_LA,
                        **kwargs,
                    )
                    self.nth_update += nth_update
                    self.delta_iota_over_run = np.concat(
                        [self.delta_iota_over_run, delta_iota_over_run]
                    )
                    self.iters_per_run = np.concat([self.iters_per_run, iters_per_run])

                    iters_used = self.get_iterations() - iters_used_start
        t1 = time.time()

        logger.info(
            f"Terminated current minimization after {t1 - t0:.2}s with Δiota = {self.rms_iota:.3e}, Δiota = {iota_tol:.3e}. Used {iters_used} / {maxIter} iterations."
        )
        self.save()

    @with_binding
    def minimize(
        self,
        iota_target: float,
        maxUpdates: int,
        target: Literal["iota", "iota_and_force", "auto"],
        **kwargs,
    ):
        if self.current_constraint:
            logger.debug("minimize with current optimization")
            self._minimize_current_constraint(
                iota_target=iota_target,
                maxUpdates=maxUpdates,
                reset_LA=self.reset_LA,
                target=target,
                **kwargs,
            )
        else:
            logger.debug("minimize with fixed iota")
            self._minimize(**kwargs)
