# Copyright (c) 2026 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
import logging
import shutil
import time
from collections.abc import Mapping
from datetime import datetime
from pathlib import Path
from typing import Literal
from abc import ABC, abstractmethod

import numpy as np
import xarray as xr

import gvec
from gvec.core.state import State
from gvec.util import CaseInsensitiveDict as cidict, compute_hfactor


DEFAULT_MINIMIZE_TOL = 1e-6  # different to default `minimize_tol` in fortran
DEFAULT_TOTALITER = 10**5  # different to default `maxIter` in fortran
AUTO_PICARD_CURRENT_IOTA_TOL_TARGET = 1e-8  # least-squares fit: limited to single precision!
DEFAULT_PICARD_CURRENT_MAXRESTARTS = 50
DIAGNOSTICS_NPOINTS = 256 + 1


class AbstractRun(ABC):
    def __init__(
        self,
        params: Mapping,
        state: Path | State | None,
        redirect_gvec_stdout: bool,
        quiet: bool,
        parameter_format: Literal["toml", "yaml"],
    ):
        """
        State of a GVEC run during a stage, e.g. a picard iteration during a current optimization run.

        Parameters
        ----------
        params : Mapping
            GVEC parameter dictionary.
        state : Path | gvec.State | None, optional
            Statefile or State to restart from. The default is None.
        redirect_gvec_stdout : bool, optional
            Whether to redirect GVEC's stdout. The default is True.
        quiet : bool, optional
            Whether to suppress all output. The default is False.

        Raises
        ------
        ValueError
            If neither params nor a state is provided.
        """

        if params:
            self.parameters = params.copy()
        elif isinstance(state, gvec.State):
            self.parameters = state.parameters.copy()
        else:
            raise ValueError(
                f"Neither valid 'params' nor a valid 'state' are provided! params: {params}, state: {state}"
            )
        self.parameters = cidict.recursive(self.parameters)

        # params: user provided parameters
        # self.parameters: global/fallback parameters
        # _state_parameters reflect the state of the Run, that is they change during each stage and/or restart
        #   after a stage they are reset to parameters. They are NOT state.parameters, as those reflect the state of a
        #   finished restart.
        if "ProjectName" not in self.parameters:
            self.parameters["ProjectName"] = "GVEC"

        self.logger = logging.getLogger("gvec.run")

        if "picard_current" not in self.parameters:
            self.parameters["picard_current"] = "off"
            self.curr_constraint = False

        match self.parameters["picard_current"]:
            case "off":
                self.curr_constraint = False
                self.parameters["picard_current"] = cidict()
            case "auto":
                self.curr_constraint = True
                self.parameters["picard_current"] = cidict()
            case cidict():
                self.curr_constraint = True

        self.filetype = parameter_format

        project_dir = Path(f"{self.parameters['ProjectName']}_gvec_stages")
        if project_dir.exists():
            self.logger.debug(f"Removing existing run directory {project_dir}")
            shutil.rmtree(project_dir)
        project_dir.mkdir()
        self.project_dir = project_dir

        self.nth_stage = 0
        self.rundir: Path = None
        self.quiet = quiet
        self.diagnostics_stage: xr.Dataset = None
        self.diagnostics_minimizer = None
        self.GVEC_iter_used = 0
        self.redirect_gvec_stdout = redirect_gvec_stdout

        self.totaliter = self.parameters.get("totaliter", DEFAULT_TOTALITER)
        if "maxIter" not in self.parameters:
            self.parameters["maxIter"] = self.totaliter
        if self.parameters["maxIter"] <= 0:
            raise ValueError(
                f"Value 'maxIter' must be greater than 0 (current value is {self.parameters['maxIter']})."
            )

        if "minimize_tol" not in self.parameters:
            self.parameters["minimize_tol"] = DEFAULT_MINIMIZE_TOL

        if "I_tor" in self.parameters:
            self.curr_constraint = True
        else:
            self.curr_constraint = False

        if "stages" in self.parameters:
            self.stages = self.parameters["stages"]
            self._has_initial_stage = False
            del self.parameters["stages"]
        else:
            # initial stage
            self._has_initial_stage = True
            self.stages = [cidict()]
            self.stages[0]["maxIter"] = 0
            self.stages[0]["init_LA"] = self.parameters.get("init_LA", True)
            self.parameters["init_LA"] = False
            # move boundary perturbation to initial stage
            if self.parameters.get("boundary_perturb", False):
                self.stages[0]["boundary_perturb"] = True
                self.parameters["boundary_perturb"] = False
            # move init from vmec to initial stage
            if self.parameters.get("whichInitEquilibrium", 0) == 1:
                self.stages[0]["whichInitEquilibrium"] = 1
                self.parameters["whichInitEquilibrium"] = 0

            # ensure at least one stage beyond initial stage
            self.stages += [cidict()]

        parameters_stages = self.parameters.copy()
        parameters_stages["stages"] = self.stages
        parameters_stages_name = (
            f"parameter_{self.parameters['ProjectName']}.stages.{self.filetype}"
        )
        self.logger.info(f"... generated {len(self.stages)} stages.")

        # load restart state
        if isinstance(state, gvec.State):
            self.logger.info(f"reading restart state from {state.statefile}")
            self.state = state

            if self.parameters.get("which_hmap", 1) != self.state.parameters.get(
                "which_hmap", 1
            ):
                self.logger.warning(
                    f"restarting with hmap={self.parameters.get('which_hmap', 1)} from hmap={self.state.parameters.get('which_hmap', 1)}."
                )

            # set initial iota profile
            if "iota" not in self.parameters:
                if "I_tor" in self.state.parameters:
                    ev = self.state.evaluate("iota", rho=np.linspace(0, 1, 51))
                    self.parameters["iota"] = {
                        "type": "interpolation",
                        "vals": ev.iota.data,
                        "rho2": (ev.rho**2).data,
                    }
                else:
                    self.parameters["iota"] = self.state.parameters["iota"]

            # compute boundary perturbation relative to restart state
            base, perturbation = gvec.util.compute_boundary_perturbation(
                self.state.parameters, self.parameters
            )
            if perturbation:
                self.parameters |= base
                if "boundary_perturb_type" not in self.parameters:
                    self.parameters["boundary_perturb_type"] = "cosm"
                self.stages[0] |= perturbation
                self.stages[0]["boundary_perturb"] = True

        # load I_tor profile (and set initial iota if not provided)
        if self.curr_constraint:
            self._set_I_tor(self.parameters)
            self.iota_rms = None
            if "iota" not in self.parameters:
                self.parameters["iota"] = cidict({"type": "polynomial", "coefs": [0.0]})
        else:
            self.I_tor = None

        # account for the change in relative path of restart, hmap, etc. files
        for key, value in self.parameters.items():
            if key.lower() in [
                "vmecwoutfile",
                "boundary_filename",
                "hmap_ncfile",
            ]:
                self.parameters[key] = Path(value).resolve()
                try:
                    if self.parameters[key].exists():
                        Path(self.parameters[key].name).symlink_to(self.parameters[key])
                    else:
                        raise FileNotFoundError(
                            f"Could not find {key} at {self.parameters[key]}"
                        )
                except FileExistsError:
                    continue

        self._state_parameters = self.parameters.copy()

        # load restart state (2nd part)
        if not isinstance(state, gvec.State):  # Path, str or None
            if state is not None:
                self.logger.warning(
                    "restarting from statefile without associated parameterfile"
                )
            initial_param_file = self.project_dir / "parameters_initial.ini"
            gvec.util.write_parameters(self._state_parameters, initial_param_file)
            self.state = gvec.State(initial_param_file, state)

        self.logger.info(f"... writing new parameters to '{parameters_stages_name}'")
        gvec.util.write_parameters(
            parameters_stages,
            project_dir / parameters_stages_name,
        )

    def _set_I_tor(self, params: Mapping):
        """Evaluate and set the target toroidal current profile at linearly spaced positions in rho.

        Raises
        ------
        gvec.errors.InvalidParameterError
            If an unknown profile type is provided.
        """
        match params["I_tor"].get("type", "polynomial"):
            case "polynomial":
                coefs = np.array(params["I_tor"]["coefs"])
                coefs *= params["I_tor"].get("scale", 1.0)
                self.I_tor = np.polynomial.Polynomial(coefs)
                if abs(coefs[0]) > 1e-8:
                    raise gvec.errors.InvalidParameterError(
                        f"Toroidal current profile not zero at magnetic axis!  I_tor(rho=0): {coefs[0]}"
                    )
            case "bspline":
                from scipy.interpolate import BSpline

                coefs = np.array(params["I_tor"]["coefs"], dtype=float)
                coefs *= params["I_tor"].get("scale", 1.0)
                knots = np.array(params["I_tor"]["knots"], dtype=float)
                deg = np.sum(knots == knots[0]) - 1
                self.I_tor = BSpline(knots, coefs, deg)
                if abs(self.I_tor(0.0)) > 1e-8:
                    raise gvec.errors.InvalidParameterError(
                        f"Toroidal current profile not zero at magnetic axis! I_tor(rho=0): {self.I_tor(0.0)}"
                    )
            case "interpolation":
                from scipy.interpolate import make_interp_spline

                y_vals = np.array(params["I_tor"]["vals"], dtype=float)
                rho2_vals = np.array(params["I_tor"]["rho2"], dtype=float)
                if min(np.sqrt(rho2_vals)) > 1e-4:
                    rho2_vals = np.append([0], rho2_vals)
                    y_vals = np.append([0], y_vals)
                self.I_tor = make_interp_spline(rho2_vals, y_vals)
                if abs(self.I_tor(0.0)) > 1e-8:
                    raise gvec.errors.InvalidParameterError(
                        f"Toroidal current profile not zero at magnetic axis! I_tor(rho=0): {self.I_tor(0.0)}"
                    )
            case _:
                raise gvec.errors.InvalidParameterError(
                    f"Unknown Itor type: {params['I_tor']['type']}"
                )

    def _reset_params_to_original(self):
        """
        Reset the parameters to the original values. Except for `iota` and `maxIter`, which is limited by `totaliter`.
        """
        params = self.parameters.copy()
        if "stages" in params:
            del params["stages"]
        params["iota"] = self._state_parameters["iota"]
        params["maxIter"] = min(
            self.totaliter - self.GVEC_iter_used,
            self.parameters.get("maxIter", self.totaliter),
        )
        self._state_parameters = params

    def _set_params_for_stage(self, stage: Mapping):
        """
        Set the run parameters to the values specified in the stage.

        Parameters
        ----------
        stage : Mapping
            Dictionary specifying which parameters are to be changed from the original parameter set.
        """
        stage = cidict(stage)
        set_I_tor = False

        for key, value in stage.lower_items():
            if key == "maxiter":
                self._state_parameters[key] = min(self.totaliter - self.GVEC_iter_used, value)
            if key == "picard_current":
                match value:
                    case "auto":
                        self._state_parameters[key] = cidict()
                        self.curr_constraint = True
                        continue
                    case "off":
                        self.curr_constraint = False
                        self._state_parameters[key] = cidict()
                        continue
                    case Mapping():
                        self.curr_constraint = True

            if key in ["iota", "pres", "sgrid", "i_tor"]:
                if key not in self._state_parameters:
                    self._state_parameters[key] = cidict()
                if key == "i_tor":
                    set_I_tor = True
                for subkey, subvalue in value.items():
                    self._state_parameters[key][subkey] = subvalue
            if key in self._state_parameters and isinstance(value, Mapping):
                for subkey, subvalue in value.items():
                    self._state_parameters[key][subkey] = subvalue
            else:
                self._state_parameters[key] = value

        if set_I_tor:
            self._set_I_tor(self._state_parameters)

    def run(
        self,
        keep_intermediates: Literal["all", "stages"] | None = None,
    ):
        """Sequentially run the stages of the Run object.

        Parameters
        ----------
        keep_intermediates : Literal["all", "stages"] | None, optional
            Whether to keep intermediate results of GVEC. With `"all"`, all intermediate results are kept. With `"stages"`,
            only the final restarts from each stage are kept. With `None`, all intermediate results are deleted. The default is None.

        Returns
        -------
        tuple:
            The final state (gvec.State) and the diagnostics dataset (xr.DataSet)

        Raises
        ------
        gvec.errors.InvalidParameterError
            If keep_intermediates is not None, "stages" or "all"
        gvec.errors.InvalidParameterError
            If stages are set when 'picard_current="auto"'
        gvec.errors.MissingParameterError
            If 'iota_tol' is not specified when 'I_tor' is provided.
        gvec.errors.InvalidParameterError
            If 'picard_current.target' is not properly specified.
        """
        if keep_intermediates and keep_intermediates not in ["all", "stages"]:
            raise gvec.errors.InvalidParameterError(
                f"""'keep_intermediates' has to be either None, "stages" or "all" but is {keep_intermediates}"""
            )
        if keep_intermediates == "all":
            self.save_intermediate = True
        else:
            self.save_intermediate = False

        start_time = time.time()
        for s, stage in enumerate(self.stages):
            if self.GVEC_iter_used >= self.totaliter:
                self.logger.warning("Maximum number of GVEC iterations reached!")
                break
            self.logger.debug(f"Entering stage {s}.")
            self.nth_stage = s

            self._reset_params_to_original()
            self._set_params_for_stage(stage)
            self._log_hfactors()

            # run the stage
            self._minimize_stage()

        self.logger.info("Done.")
        end_time = time.time()

        final_iota_str = ""
        if self.curr_constraint:
            final_iota_str += f"\n and rms Δiota = {self.state.rms_iota:.2e}"

        final_message = (
            f"GVEC finished after {end_time - start_time:5.1f} seconds",
            f" using {self.GVEC_iter_used} iterations (totalIter = {self.totaliter})",
            f" with |force| = {self.max_force:.2e} (minimize_tol = {self._state_parameters['minimize_tol']:.2e})",
            final_iota_str,
        )
        if self.quiet:
            self.logger.info("".join(final_message))
        else:
            print(*final_message)
        final_statefile = Path(self._state_parameters["ProjectName"] + "_State_final.dat")
        final_parameter_file = Path(
            "parameter_" + self._state_parameters["ProjectName"] + "_final.ini"
        )

        shutil.copy(self.state.statefile, final_statefile)
        parameters_final = self.state._parameter_dict

        parameters_final["maxIter"] = -1
        for key in parameters_final:
            if key.lower() in [
                "vmecwoutfile",
                "boundary_filename",
                "hmap_ncfile",
            ]:
                parameters_final[key] = Path(self.parameters[key]).name
        parameters_final = gvec.util.flatten_parameters(parameters_final)
        gvec.util.write_parameter_file_ini(
            parameters=parameters_final,
            path=final_parameter_file,
            header=f"!Auto-generated with `pygvec run` (stage {self.nth_stage})\n"
            f"!Created at {datetime.now().isoformat()}\n"
            f"!pyGVEC v{gvec.__version__}\n",
        )
        self.state = gvec.State(final_parameter_file.absolute(), final_statefile.absolute())
        if keep_intermediates is None:
            shutil.rmtree(self.project_dir)
            self.project_dir = None

        diagnostics = (self.diagnostics_stage, self.diagnostics_minimizer)
        return (self.state, diagnostics)

    def _minimize_stage(self):
        self._state_parameters["maxIter"] = min(
            self.totaliter - self.GVEC_iter_used,
            self._state_parameters["maxIter"],
        )
        self._print_progress()

        self._run_single_minimization()
        if self.curr_constraint:
            self._state_parameters["iota"] = self.state._parameter_dict["iota"]
        if (
            self.max_force > self._state_parameters["minimize_tol"]
            and self._state_parameters["maxIter"] > 0
        ):
            self.logger.warning(
                f"Force tolerance was not reached! |force|: {self.max_force:.2e}, minimize_tol: {self._state_parameters['minimize_tol']:.2e}"
            )

    def _print_progress(self):
        """
        Evaluate and print a progress string for the current stage and run.
        """

        progressstr = "|"
        progressstr += "=|" * self.nth_stage
        progressstr += ">|" + ".|" * (len(self.stages) - self.nth_stage - 1)
        start_str = "GVEC"
        state_str = ""
        if len(self.stages) > 1:
            state_str += f" - completed {self.nth_stage}/{len(self.stages)} stages"
        progressstr = start_str + state_str + progressstr
        if self.logger.isEnabledFor(logging.INFO):
            self.logger.info(progressstr)
        elif not self.quiet:
            print(progressstr, end="\r")

    def _log_hfactors(self):
        try:
            M1, N1 = self._state_parameters["X1_mn_max"]
            M2, N2 = self._state_parameters["X2_mn_max"]
            Ml, Nl = self._state_parameters["LA_mn_max"]
            K = self._state_parameters["sgrid"]["nelems"]
            if (M1 == M2 == Ml) and (N1 == N2 == Nl):
                hM, hN, hK = compute_hfactor(M1, N1, K)
                self.logger.info(f"h-factors: pol={hM:.2e}, tor={hN:.2e}, rad={hK:.2e}")
            else:
                hM1, hN1, hK1 = compute_hfactor(M1, N1, K)
                hM2, hN2, _ = compute_hfactor(M2, N2, K)
                hMl, hNl, _ = compute_hfactor(Ml, Nl, K)
                self.logger.info(
                    f"h-factors: pol(X1)={hM1:.2e}, pol(X2)={hM2:.2e}, pol(LA)={hMl:.2e}, tor(X1)={hN1:.2e}, tor(X2)={hN2:.2e}, tor(LA)={hNl:.2e}, rad={hK1:.2e}"
                )
        except (KeyError, TypeError, ValueError):
            self.logger.debug(
                "Could not compute h-factors due to missing or invalid parameters."
            )

    @abstractmethod
    def _run_single_minimization(self):
        """Runs a single minimization within a stage."""
