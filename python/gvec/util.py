# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""GVEC utility module

This module is part of the gvec python package, but also used directly in the tests.
"""

import contextlib
import copy
import os
import re
import shutil
from collections.abc import Mapping, MutableMapping, Iterable
from pathlib import Path
from typing import Literal
from copy import deepcopy
import logging

import numpy as np
from numpy.typing import ArrayLike

try:
    from scipy.interpolate import BSpline
except ImportError:
    BSpline = None

logger = logging.getLogger(__name__)


def get_compile_options() -> str:
    try:
        from gvec import _compile_options as opts
    except ImportError:
        return "UNKNOWN BUILD"
    config = opts.CMAKE_BUILD_TYPE
    if opts.USE_OPENMP:
        config += ", OpenMP"
    if opts.USE_MPI:
        config += ", MPI"
    if opts.GVEC_FIX_HMAP:
        config += f", {opts.GVEC_FIX_HMAP}"
    config += f", {opts.CMAKE_Fortran_COMPILER.name}"
    if len(opts.CMAKE_HOSTNAME) > 0:
        config += f" on {opts.CMAKE_HOSTNAME}"
    return config


def version_info() -> str:
    import platform
    from gvec._version import __version__

    return f"pyGVEC v{__version__} ({get_compile_options()}, python {platform.python_version()}) from {Path(__file__).parent}"


@contextlib.contextmanager
def chdir(target: Path | str):
    """
    Contextmanager to change the current working directory.

    Using a context has the benefit of automatically changing back to the original directory when the context is exited, even if an exception is raised.
    """
    target = Path(target)
    source = Path.cwd()

    try:
        os.chdir(target)
        yield
    finally:
        os.chdir(source)


class CaseInsensitiveDict(MutableMapping):
    # Adapted from requests.structures.CaseInsensitiveDict
    # See: https://github.com/psf/requests/blob/main/src/requests/structures.py
    # Original license: Apache License 2.0
    """A dictionary-like Mutable Mapping where string keys are case-insensitive.

    Implements all methods and operations of
    ``MutableMapping`` as well as dict's ``copy``. Also
    provides ``lower_items`` and ``lower_keys``.

    Keys that are not strings will be stored as-is.
    The structure remembers the case of the last used key, and
    ``iter(instance)``, ``keys()``, ``items()``, ``iterkeys()``, and ``iteritems()``
    will contain case-sensitive keys. However, querying and contains
    testing is case insensitive:

        cid = CaseInsensitiveDict()
        cid['param'] = 'value'
        cid['Param'] == 'value'  # True

    If the constructor, ``.update``, or equality comparison
    operations are given keys that have equal ``.lower()``s, a ValueError is raised.
    """

    def __init__(self, data=(), /, **kwargs):
        self._data = {}
        self.update(data, **kwargs)

    @staticmethod
    def _idx(key):
        return key.lower() if isinstance(key, str) else key

    def __setitem__(self, key, value):
        # Use the lowercased key for lookups, but remember the last key alongside the value.
        self._data[self._idx(key)] = (key, value)

    def __getitem__(self, key):
        return self._data[self._idx(key)][1]

    def __delitem__(self, key):
        del self._data[self._idx(key)]

    def update(self, data=(), /, **kwargs):
        updates = {}
        updates.update(data, **kwargs)
        idxs = {self._idx(key) for key in updates}
        if len(idxs) != len(updates):
            raise ValueError("Duplicate keys passed to CaseInsensitiveDict.update")
        for key, value in updates.items():
            self[key] = value

    def __iter__(self):
        return (key for key, value in self._data.values())

    def lower_keys(self):
        return (idx for idx in self._data.keys())

    def lower_items(self):
        return ((idx, value) for idx, (key, value) in self._data.items())

    def __len__(self):
        return len(self._data)

    def __eq__(self, other):
        if isinstance(other, Mapping):
            other = CaseInsensitiveDict(other)
        else:
            return NotImplemented
        # Compare insensitively
        return dict(self.lower_items()) == dict(other.lower_items())

    def serialize(self):
        """Recursively serialize this object, converting Mappings to dicts and Iterables to lists."""

        def _serialize(value):
            if isinstance(value, Mapping):
                return {k: _serialize(v) for k, v in value.items()}
            elif isinstance(value, Iterable) and not isinstance(value, str):
                return [_serialize(v) for v in value]
            else:
                return value

        return _serialize(self)

    def __repr__(self):
        return f"{self.__class__.__name__}{dict(self.items())}"

    def copy(self):
        """Return a deep copy."""
        return deepcopy(self)

    def __or__(self, other):
        """Union/Merge operator 'a | b' (without modifying self)."""
        if not isinstance(other, Mapping):
            return NotImplemented
        result = CaseInsensitiveDict(self)
        result.update(other)
        return result

    def __ior__(self, other):
        """In-place union/merge operator 'a |= b' (modifies self)."""
        if not isinstance(other, Mapping):
            return NotImplemented
        self.update(other)
        return self


def adapt_parameter_file(source: str | Path, target: str | Path, **kwargs):
    """
    Copy the `source` file to the `target` file and replace the parameters according to `kwargs`.

    Args:
        source (str or Path): The path to the source parameter file.
        target (str or Path): The path to the target parameter file.
        **kwargs: Keyword arguments representing the parameters to be replaced.
                  if the value of the key is "!", the line with the keyword is uncommented, if possible

    Raises:
        AssertionError: If the number of occurrences for any parameter is not exactly 1.

    Notes:
        - If no parameters are provided in `kwargs`, the function simply copies the `source` file to the `target` file.
        - The function replaces the parameters in the format `key = value`, where value is either a sequence of characters containing
          no whitespace or a single pair of parentheses with any content. The value from `kwargs` is inserted using the standard python
          string conversion. There may be a comment, starting with `!`, after the value.
        - If a parameter already exists in the `source` file, its value is replaced with the corresponding value from `kwargs`.
        - If a parameter does not exist in the `source` file, it is added to the `target` file.
        - If the value of the key starts with "!", the line with the keyword is just uncommented.  (i.e. "!key=2.5" -> "key=2.5")
          If no line with the keyword is found, the key is added with the value, excluding the leading "!"  (i.e. value is "!0.5" -> "key=0.5" is added)

    Example:
        `>>> adapt_parameter_file('/path/to/source.ini', '/path/to/target.ini', param1=1.2, param2="(1, 2, 3)")`
    """
    if not len(kwargs.keys()):
        shutil.copy2(source, target)
        return

    for key, value in kwargs.items():
        if isinstance(value, Mapping) or isinstance(value, str):
            pass
        elif isinstance(value, bool):
            kwargs[key] = "T" if value else "F"
        elif isinstance(value, Iterable):
            kwargs[key] = f"(/{', '.join(map(str, value))}/)"
        else:
            kwargs[key] = str(value)
    kwargs = {key.lower(): value for key, value in kwargs.items()}

    # initialize occurrences counters for all parameters to be set
    occurrences = {}
    for key in kwargs:
        if isinstance(kwargs[key], Mapping):
            for m, n in kwargs[key]:
                occurrences[key, m, n] = 0
        else:
            occurrences[key] = 0

    with open(source, "r") as source_file:
        source_file = source_file.readlines()
    with open(target, "w") as target_file:
        for line in source_file:
            if m := re.match(
                r"\s*([^!=\s\(]+)\s*\(\s*([-\d]+);\s*([-\d]+)\)\s*=\s*([-+\d\.Ee]+)",
                line,
            ):
                key, *mn, value = m.groups()
                if key.lower() in kwargs:
                    if (int(mn[0]), int(mn[1])) in kwargs[key.lower()]:
                        line = f"{key}({mn[0]};{mn[1]}) = {kwargs[key.lower()][(int(mn[0]), int(mn[1]))]}\n"
                        occurrences[key.lower(), int(mn[0]), int(mn[1])] += 1
            elif m := re.match(
                r"([\s!]*)("
                + "|".join(
                    [
                        key.lower()
                        for key, value in kwargs.items()
                        if not isinstance(value, Mapping)
                    ]
                )
                + r")(\s*=\s*)(\([^\)]*\)|[^!\s]*)(.*)",
                line,
                re.IGNORECASE,
            ):
                prefix, key, sep, value, suffix = m.groups()
                if "!" in prefix:  # found commented keyword
                    if str(kwargs[key.lower()])[0] == "!":  # only uncomment keyword
                        line = f"{key}{sep}{value}{suffix}\n"
                        occurrences[key.lower()] += 1
                else:  # found uncommented keywords
                    if str(kwargs[key.lower()])[0] != "!":  # use new keyword
                        line = f"{prefix}{key}{sep}{kwargs[key.lower()]}{suffix}\n"
                        occurrences[key.lower()] += 1
                    else:  # use the existing keyword,value pair with a comment
                        line = (
                            f"{prefix}{key}{sep}{value} !!WAS ALREADY UNCOMMENTED!! {suffix}\n"
                        )
                        occurrences[key.lower()] += 1
            target_file.write(line)
        # add key,value pair if not existing in parameterfile.
        for key, o in occurrences.items():
            if o == 0:
                if isinstance(key, tuple):
                    key, m, n = key
                    if str(kwargs[key][m, n]) != "!":
                        target_file.write(f"\n{key}({m};{n}) = {kwargs[key][m, n]}")
                        occurrences[key, m, n] += 1
                else:
                    if str(kwargs[key]) == "!":
                        continue  # ignore 'uncomment' value if key is not found
                    elif str(kwargs[key])[0] == "!":
                        # use default value '!default' if key is not found
                        target_file.write(f"\n{key} = {kwargs[key][1:]}")
                    else:
                        # add parameter at the end if key is not found
                        target_file.write(f"\n{key} = {kwargs[key]}")
                    occurrences[key] += 1
    assert all([o == 1 for o in occurrences.values()]), (
        f"bad number of occurrences in adapt_parameter_file: {occurrences}"
    )


def write_parameter_file_ini(
    parameters: Mapping, path: str | Path = "parameter.ini", header: str = ""
):
    """
    Write the parameters to the specified parameter file in GVEC-ini format.

    Args:
        parameters: A mapping containing the parameters to be written to the parameter file.
        path: The path to the parameter file.
    """
    parameters = parameters.copy()
    for key, value in parameters.items():
        if isinstance(value, Mapping) or isinstance(value, str):
            pass
        elif isinstance(value, bool):
            parameters[key] = "T" if value else "F"
        elif isinstance(value, Iterable):
            parameters[key] = f"(/{', '.join(map(str, value))}/)"
        else:
            parameters[key] = str(value)

    with open(path, "w") as file:
        file.write(header)
        for key, value in parameters.items():
            if isinstance(value, Mapping):
                for (m, n), val in value.items():
                    file.write(f"{key}({m};{n}) = {val}\n")
            else:
                file.write(f"{key} = {value}\n")


def read_parameter_file_ini(path: str | Path) -> CaseInsensitiveDict:
    """
    Read the parameters from the specified parameter file in GVEC-ini format.

    Args:
        path (str | Path): The path to the parameter file.

    Returns:
        CaseInsensitiveDict: A mapping (with case insensitive keys) containing the parameters from the parameter file.

    Example:
    >>> read_parameter_file_ini('/path/to/parameter.ini')
    {'param1': 1.2, 'param2': (1, 2, 3), 'param3': {(-1, 0): 0.5, (0, 0): 1.0}}
    """
    INT = r"[-+]?\d+"
    FLOAT = r"[-+]?\d*\.?\d*(?:[eE][-+]?\d+)?"
    STR = r"\S+"
    KEY = r"\w+"

    def convert(value: str):
        if "," in value:
            return tuple(convert(v) for v in value.split(","))
        if re.fullmatch(INT, value):
            return int(value)
        if re.fullmatch(FLOAT, value):
            return float(value)
        if value.upper() == "T":
            return True
        if value.upper() == "F":
            return False
        if re.fullmatch(STR, value):
            return value
        raise ValueError(f"Cannot parse value '{value}' in parameter file {path}")

    # follow the implementation in src/globals/readintools.f90:FillStrings
    parameters = CaseInsensitiveDict()
    with open(path, "r") as file:
        # read lines and preprocess them
        lines = []
        for line in file:
            # remove comments `!` and `#`
            line = re.split(r"[!#]", line)[0]
            # remove array brackets `(/` and `/)`
            line = re.sub(r"\(\/", "", line)
            line = re.sub(r"\/\)", "", line)
            # remove whitespace
            line = re.sub(r"\s+", "", line).strip()
            # skip empty lines
            if len(line) == 0:
                continue
            # combine lines that end with a `&`
            if lines and lines[-1].endswith("&"):
                lines[-1] = lines[-1][:-1] + line
            else:
                lines.append(line)

        # parse the lines
        for line in lines:
            # match parameter in the form `key(m;n) = value` with m,n integers
            if ln := re.fullmatch(rf"({KEY})\(({INT});({INT})\)=(.+)", line):
                key, m, n, value = ln.groups()
                m, n = int(m), int(n)

                if key in parameters and not isinstance(parameters[key], MutableMapping):
                    raise TypeError(
                        f"Trying to set indices for parameter '{key}' in {path}, but it is already set to a non-mapping value: {parameters[key]}"
                    )
                if key not in parameters:
                    parameters[key] = {}
                if (m, n) in parameters[key]:
                    raise IndexError(
                        f"Duplicate indices ({m}, {n}) for parameter '{key}' in {path}"
                    )
                parameters[key][m, n] = convert(value)
            # match parameter in the form `key = value`
            elif "=" in line:
                key, value = line.split("=", 1)
                if key in parameters:
                    raise IndexError(f"Duplicate parameter '{key}' in {path}")
                if not re.fullmatch(KEY, key):
                    raise ValueError(f"Invalid key '{key}' in parameter file {path}")
                parameters[key] = convert(value)
    return parameters


def check_boundary_direction(parameters: Mapping) -> bool:
    """Determine whether the boundary is described by right-handed logical coordinates (θ,ζ).

    GVEC requires a right-handed logical coordinate system (ρ,θ,ζ).
    The logical coordinate system of the poloidal plane, (ρ,θ) is also required to be right-handed,
    which requires the poloidal angle to increase in the counter-clockwise direction.
    As a consequence the toroidal angle has to increase in the clockwise direction when viewed from above.
    This is ensured in the definition of the h-maps.

    Returns:
        bool: True if (ρ,θ) is right-handed / θ increases counter-clockwise, False otherwise.
    """
    return signed_cross_sectional_area(parameters, 0.0) > 0


def signed_cross_sectional_area(
    parameters: Mapping, zeta: float, resolution: int = 1000
) -> float:
    t = np.linspace(0, 2 * np.pi, resolution, endpoint=False)
    x1 = np.zeros_like(t)
    dx1dt = np.zeros_like(t)
    x2 = np.zeros_like(t)
    dx2dt = np.zeros_like(t)
    nfp = parameters.get("nfp", 1)
    for (m, n), value in parameters.get("X1_b_cos", {}).items():
        x1 += value * np.cos(m * t - n * nfp * zeta)
        dx1dt -= value * m * np.sin(m * t - n * nfp * zeta)
    for (m, n), value in parameters.get("X1_b_sin", {}).items():
        x1 += value * np.sin(m * t - n * nfp * zeta)
        dx1dt += value * m * np.cos(m * t - n * nfp * zeta)
    for (m, n), value in parameters.get("X2_b_cos", {}).items():
        x2 += value * np.cos(m * t - n * nfp * zeta)
        dx2dt -= value * m * np.sin(m * t - n * nfp * zeta)
    for (m, n), value in parameters.get("X2_b_sin", {}).items():
        x2 += value * np.sin(m * t - n * nfp * zeta)
        dx2dt += value * m * np.cos(m * t - n * nfp * zeta)
    dA = x1 * dx2dt - x2 * dx1dt
    return np.sum(dA)


def effective_minor_radius(
    parameters: Mapping,
    resolution: tuple[int, int] = (1000, 100),
):
    nfp = parameters.get("nfp", 1)
    areas = np.zeros(resolution[1])
    for z, zeta in enumerate(np.linspace(0, 2 * np.pi / nfp, resolution[1], endpoint=False)):
        areas[z] = abs(signed_cross_sectional_area(parameters, zeta, resolution=resolution[0]))
    return np.sqrt(np.mean(areas) / np.pi)


def evaluate_boundary(
    theta: np.ndarray, zeta: np.ndarray, parameters: Mapping
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the boundary at the given (theta, zeta) points.

    Args:
        theta (1D np.ndarray): The poloidal angles at which to evaluate the boundary.
        zeta (1D np.ndarray): The toroidal angles at which to evaluate the boundary.
        parameters (Mapping): The parameters defining the boundary.

    Returns:
        tuple[2D np.ndarray, 2D np.ndarray]: The (X^1, X^2) coordinates of the boundary at the given (theta, zeta) points.
    """
    theta = np.asarray(theta)
    zeta = np.asarray(zeta)
    if theta.ndim != 1 or zeta.ndim != 1:
        raise ValueError("theta and zeta must be 1D arrays")
    nfp = parameters.get("nfp", 1)
    theta, zeta = np.meshgrid(theta, zeta, indexing="ij")
    x1 = np.zeros_like(theta)
    x2 = np.zeros_like(theta)
    for (m, n), value in parameters.get("X1_b_cos", {}).items():
        x1 += value * np.cos(m * theta - n * nfp * zeta)
    for (m, n), value in parameters.get("X1_b_sin", {}).items():
        x1 += value * np.sin(m * theta - n * nfp * zeta)
    for (m, n), value in parameters.get("X2_b_cos", {}).items():
        x2 += value * np.cos(m * theta - n * nfp * zeta)
    for (m, n), value in parameters.get("X2_b_sin", {}).items():
        x2 += value * np.sin(m * theta - n * nfp * zeta)
    return x1, x2


def evaluate_axis(zeta: np.ndarray, parameters: Mapping) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the magnetic axis at the given zeta points.

    Args:
        zeta (1D np.ndarray): The toroidal angles at which to evaluate the axis.
        parameters (Mapping): The parameters defining the axis.

    Returns:
        tuple[1D np.ndarray, 1D np.ndarray]: The (X^1, X^2) coordinates of the axis at the given zeta points.
    """
    zeta = np.asarray(zeta)
    if zeta.ndim != 1:
        raise ValueError("zeta must be a 1D array")
    nfp = parameters.get("nfp", 1)
    x1 = np.zeros_like(zeta)
    x2 = np.zeros_like(zeta)
    for (m, n), value in parameters.get("X1_a_cos", {}).items():
        if m != 0:
            raise ValueError("Axis X1_a_cos should only have m=0 modes")
        x1 += value * np.cos(-n * nfp * zeta)
    for (m, n), value in parameters.get("X1_a_sin", {}).items():
        if m != 0:
            raise ValueError("Axis X1_a_sin should only have m=0 modes")
        x1 += value * np.sin(-n * nfp * zeta)
    for (m, n), value in parameters.get("X2_a_cos", {}).items():
        if m != 0:
            raise ValueError("Axis X2_a_cos should only have m=0 modes")
        x2 += value * np.cos(-n * nfp * zeta)
    for (m, n), value in parameters.get("X2_a_sin", {}).items():
        if m != 0:
            raise ValueError("Axis X2_a_sin should only have m=0 modes")
        x2 += value * np.sin(-n * nfp * zeta)
    return x1, x2


def compute_boundary_perturbation(
    base_parameters: Mapping, perturbed_parameters: Mapping
) -> tuple[CaseInsensitiveDict, CaseInsensitiveDict]:
    """Computes the difference between the perturbed and base boundary parameters as a perturbation."""
    new_base = CaseInsensitiveDict()
    new_perturbed = CaseInsensitiveDict()
    for i in [1, 2]:
        for sincos in ["sin", "cos"]:
            perturbed = {}
            base = {}
            # set boundary modes to values from restart
            for (m, n), v in base_parameters.get(f"{i}_b_{sincos}", {}).items():
                base[m, n] = v
            # set boundary perturbation to difference between current and restart
            for (m, n), v in perturbed_parameters.get(f"X{i}_b_{sincos}", {}).items():
                v = v - base_parameters.get(f"X{i}_b_{sincos}", {}).get((m, n), 0)
                if v != 0.0:
                    perturbed[m, n] = v
            if base or perturbed:
                new_base[f"X{i}_b_{sincos}"] = base
                new_perturbed[f"X{i}pert_b_{sincos}"] = perturbed
    return new_base, new_perturbed


def flip_boundary_theta(parameters: MutableMapping) -> MutableMapping:
    """Flip the boundary parameters in the poloidal direction. θ → -θ."""
    output_params = copy.deepcopy(parameters)
    for var in ["X1_b", "X2_b"]:
        if f"{var}_cos" in parameters:
            output_params[f"{var}_cos"] = {}
            for (m, n), value in parameters[f"{var}_cos"].items():
                if m == 0:
                    output_params[f"{var}_cos"][m, n] = value
                else:
                    output_params[f"{var}_cos"][m, -n] = value
        if f"{var}_sin" in parameters:
            output_params[f"{var}_sin"] = {}
            for (m, n), value in parameters[f"{var}_sin"].items():
                if m == 0:
                    output_params[f"{var}_sin"][m, n] = value
                else:
                    output_params[f"{var}_sin"][m, -n] = -value
    return output_params


def flip_boundary_zeta(parameters: MutableMapping) -> MutableMapping:
    output_params = copy.deepcopy(parameters)
    for var in ["X1_b", "X2_b", "X1_a", "X2_a"]:
        if f"{var}_cos" in parameters:
            output_params[f"{var}_cos"] = {}
            for (m, n), value in parameters[f"{var}_cos"].items():
                if m == 0:
                    output_params[f"{var}_cos"][m, n] = value
                else:
                    output_params[f"{var}_cos"][m, -n] = value
        if f"{var}_sin" in parameters:
            output_params[f"{var}_sin"] = {}
            for (m, n), value in parameters[f"{var}_sin"].items():
                if m == 0:
                    output_params[f"{var}_sin"][m, n] = -value
                else:
                    output_params[f"{var}_sin"][m, -n] = value
    return output_params


def shift_boundary_theta_pi(parameters: MutableMapping) -> MutableMapping:
    """
    Shift the theta origin of the boundary by pi.

    cos(m(θ+π)-nζ) = (-1)^m cos(mθ-nζ)
    sin(m(θ+π)-nζ) = (-1)^m sin(mθ-nζ)
    """
    parameters = copy.deepcopy(parameters)
    for var in ["X1", "X2"]:
        for sc in ["cos", "sin"]:
            key = f"{var}_b_{sc}"
            if key in parameters:
                for (m, n), value in list(parameters[key].items()):
                    if m % 2 == 1:
                        parameters[key][m, n] = -value
    return parameters


def flip_parameters_theta(parameters: MutableMapping) -> MutableMapping:
    parameters = flip_boundary_theta(parameters)

    for profile in ["iota"]:
        if profile in parameters:
            parameters[profile]["scale"] = -parameters[profile].get("scale", 1.0)

    return parameters


def flip_parameters_zeta(parameters: MutableMapping) -> MutableMapping:
    parameters = flip_boundary_zeta(parameters)

    if "phiedge" in parameters:
        parameters["phiedge"] = -parameters["phiedge"]
    for profile in ["iota", "I_tor"]:
        if profile in parameters:
            parameters[profile]["scale"] = -parameters[profile].get("scale", 1.0)

    return parameters


def parameters_from_vmec(nml: Mapping, name: str) -> CaseInsensitiveDict:
    def as_list(value) -> list:
        if isinstance(value, list):
            return value
        else:
            return [value]

    try:
        nml = nml.todict()  # f90nml.Namelist -> dict | fills '_start_index' attribute
    except AttributeError:
        pass

    M, N = nml.get("mpol", 2) - 1, nml.get("ntor", 0)
    stellsym = not nml.get("lasym", False)  # stellarator symmetry
    params = CaseInsensitiveDict(
        ProjectName=name,
        which_hmap=1,
        minimize_tol=1e-7,
        totalIter=10000,
        logIter=100,
        nfp=nml.get("nfp", 1),
        PhiEdge=nml.get("phiedge", 1.0),
        X1X2_deg=5,
        LA_deg=5,
        sgrid=dict(
            grid_type=0,
            nElems=5,
        ),
    )
    # --- profiles --- #
    if nml.get("pmass_type", "power_series") != "power_series":
        raise ValueError(
            f"VMEC pressure profile of type {nml['pmass_type']} is not supported for conversion"
        )
    if "am" in nml:
        params["pres"] = {
            "type": "polynomial",
            "coefs": as_list(nml["am"]),
            "scale": nml.get("pres_scale", 1.0),
        }
    else:
        logger.warning("No pressure profile defined.")
    if nml.get("piota_type", "power_series") != "power_series":
        raise ValueError(
            f"VMEC iota profile of type {nml['piota_type']} is not supported for conversion"
        )
    if "ai" in nml:
        params["iota"] = {
            "type": "polynomial",
            "coefs": as_list(nml["ai"]),
        }
    if nml.get("ncurr", 0) == 1:  # ncurr = 0: flux conservation | ncurr = 1: current constraint
        if nml.get("pcurr_type", "power_series") != "power_series":
            raise ValueError(
                f"VMEC current profile of type {nml['pcurr_type']} is not supported for conversion"
            )
        params["I_tor"] = {"type": "polynomial"}
        if nml["curtor"] == 0.0:
            params["I_tor"]["coefs"] = [0.0]
            params["I_tor"]["scale"] = 1.0
        else:
            coefs = [0] + [
                p / (i + 1) for i, p in enumerate(as_list(nml["ac"]))
            ]  # I'(s) -> I(s)
            params["I_tor"]["coefs"] = coefs
            params["I_tor"]["scale"] = nml["curtor"] / sum(coefs)
        params["picard_current"] = "auto"
    if "ai" not in nml and nml.get("ncurr", 0) == 0:
        logger.warning("No iota or current profile defined.")

    # --- boundary --- #
    for vmec_key, gvec_key in [
        ("rbc", "X1_b_cos"),
        ("rbs", "X1_b_sin"),
        ("zbc", "X2_b_cos"),
        ("zbs", "X2_b_sin"),
    ]:
        if vmec_key not in nml:
            continue
        values = np.array(nml[vmec_key], dtype=float)
        if "_start_index" in nml and vmec_key in nml["_start_index"]:
            n0, m0 = nml["_start_index"][vmec_key]
            M = max(M, values.shape[0] - 1 + m0)
            N = max(N, abs(n0), values.shape[1] - 1 + max(0, n0))
        else:
            raise ValueError(
                f"VMEC namelist array '{vmec_key}' has shape {values.shape} that does not match the expected shape {(M + 1, 2 * N + 1)=} and no '_start_index' is given in the namelist."
            )
        params[gvec_key] = {}
        for m in range(m0, m0 + values.shape[0]):
            for n in range(n0, n0 + values.shape[1]):
                if m == 0 and n < 0:
                    continue
                v = values[m - m0, n - n0]
                if not np.isnan(v) and v != 0.0:
                    params[gvec_key][m, n] = values[m - m0, n - n0]

    if "rbs" in nml or "zbc" in nml:
        if stellsym:
            logger.warning(
                "VMEC namelist contains 'RBS' or 'ZBC' but is supposed to be stellarator symmetric. Assuming asymmetry."
            )
        params["X1_sin_cos"] = "_sincos_"
        params["X2_sin_cos"] = "_sincos_"
        params["LA_sin_cos"] = "_sincos_"
    else:
        if not stellsym:
            logger.warning(
                "VMEC namelist does not contain 'RBS' or 'ZBC' but is supposed to be non-stellarator symmetric. Assuming symmetry."
            )
        params["X1_sin_cos"] = "_cos_"
        params["X2_sin_cos"] = "_sin_"
        params["LA_sin_cos"] = "_sin_"

    # --- axis --- #
    for vmec_key, gvec_key in [
        ("raxis_cc", "X1_a_cos"),
        ("raxis_cs", "X1_a_sin"),
        ("zaxis_cc", "X2_a_cos"),
        ("zaxis_cs", "X2_a_sin"),
    ]:
        if vmec_key not in nml:
            continue
        values = np.array(nml[vmec_key], dtype=float)
        if values.ndim == 0:
            continue
        params[gvec_key] = {(0, n): v for n, v in enumerate(values)}
    if not any(k in params for k in ["X1_a_cos", "X1_a_sin", "X2_a_cos", "X2_a_sin"]):
        params["init_average_axis"] = True

    # --- other --- #
    params.update(
        X1_mn_max=(M, N),
        X2_mn_max=(M, N),
        LA_mn_max=(M, N),
    )
    return params


def axis_from_boundary(parameters: MutableMapping) -> MutableMapping:
    parameters2 = copy.deepcopy(parameters)
    N = parameters["X1_mn_max"][1]
    parameters2["X1_a_cos"] = {(0, n): parameters["X1_b_cos"][0, n] for n in range(N + 1)}
    parameters2["X2_a_sin"] = {(0, n): parameters["X2_b_sin"][0, n] for n in range(N + 1)}
    if "X1_b_sin" in parameters:
        parameters2["X1_a_sin"] = {(0, n): parameters["X1_b_sin"][0, n] for n in range(N + 1)}
    if "X2_b_cos" in parameters:
        parameters2["X2_a_cos"] = {(0, n): parameters["X2_b_cos"][0, n] for n in range(N + 1)}
    return parameters2


def stack_parameters(parameters: Mapping) -> CaseInsensitiveDict:
    """Stack parameters into a hierarchical dictionary"""
    output = CaseInsensitiveDict()
    for key, value in parameters.items():
        if "_" not in key:
            output[key] = value
            continue
        group, name = key.split("_", 1)
        if group in ["iota", "pres", "sgrid"]:
            if group not in output:
                output[group] = CaseInsensitiveDict()
            output[group][name] = value
        else:
            output[key] = value
    return output


def flatten_parameters(parameters: Mapping) -> CaseInsensitiveDict:
    """Flatten parameters from a hierarchical dictionary"""
    output = CaseInsensitiveDict()
    for key, value in parameters.items():
        if key.lower() in ["stages", "i_tor", "picard_current", "totaliter"]:
            continue  # not supported by fortran-GVEC
        elif isinstance(value, Mapping) and not re.match(
            r"(x1|x2|la)(pert:?)?_[a|b]_(sin|cos)", key.lower()
        ):
            for subkey, subvalue in value.items():
                output[f"{key}_{subkey}"] = subvalue
        else:
            output[key] = value
    return output


def stringify_mn_parameters(parameters: Mapping) -> CaseInsensitiveDict:
    """Serialize parameters into a string"""
    output = CaseInsensitiveDict()
    for key, value in parameters.items():
        if re.match(r"(x1|x2|la)(pert:?)?_[a|b]_(sin|cos)", key.lower()):
            output[key] = {}
            for (m, n), val in value.items():
                if isinstance(val, np.number):
                    val = val.item()
                output[key][f"({m}, {n:2d})"] = val
        elif key.lower() == "stages":
            output[key] = [stringify_mn_parameters(stage) for stage in value]
        elif isinstance(value, np.number) or (
            isinstance(value, np.ndarray) and value.size == 1
        ):
            output[key] = value.item()
        elif isinstance(value, np.ndarray):
            output[key] = value.tolist()
        else:
            output[key] = value
    return output


def unstringify_mn_parameters(parameters: Mapping) -> CaseInsensitiveDict:
    """Deserialize parameters from a string"""
    output = CaseInsensitiveDict()
    for key, value in parameters.items():
        if re.match(r"(x1|x2|la)(pert:?)?_[a|b]_(sin|cos)", key.lower()):
            output[key] = CaseInsensitiveDict()
            for mn, val in value.items():
                m, n = map(int, mn.strip("()").split(","))
                output[key][(m, n)] = val
        elif key.lower() == "stages":
            output[key] = [unstringify_mn_parameters(stage) for stage in value]
        else:
            output[key] = value
    return output


def read_parameters(
    path: Path | str, format: Literal["ini", "yaml", "toml"] | None = None
) -> CaseInsensitiveDict:
    import tomlkit
    import yaml

    path = Path(path)
    # auto-detect format
    if format is None:
        format = path.suffix[1:]

    if format == "ini":
        inputs = read_parameter_file_ini(path)
        inputs = stack_parameters(inputs)
    elif format == "yaml":
        with open(path, "r") as file:
            inputs = yaml.safe_load(file)
        inputs = unstringify_mn_parameters(inputs)
    elif format == "toml":
        with open(path, "r") as file:
            inputs = tomlkit.parse(file.read()).unwrap()
        inputs = unstringify_mn_parameters(inputs)
    else:
        raise ValueError(f"Unknown parameter file format {format}")
    return inputs


def write_parameters(
    parameters: Mapping,
    path: Path | str = "parameter.ini",
    format: Literal["ini", "yaml", "toml"] | None = None,
):
    import tomlkit
    import yaml

    path = Path(path)
    # auto-detect format
    if format is None:
        format = path.suffix[1:]

    if format == "ini":
        outputs = flatten_parameters(parameters)
        write_parameter_file_ini(outputs, path)
    elif format == "yaml":
        outputs = stringify_mn_parameters(parameters)
        with open(path, "w") as file:
            yaml.safe_dump(
                outputs.serialize(), file, sort_keys=False
            )  # ToDo: specify style/flow?
    elif format == "toml":
        outputs = stringify_mn_parameters(parameters)
        with open(path, "w") as file:
            file.write(
                tomlkit.dumps(outputs.serialize())
            )  # ToDo: nicer output using document API
    else:
        raise ValueError(f"Unknown parameter file format {format}")


def bspl2gvec(
    name: Literal["iota", "pres"],
    bspl: BSpline = None,
    knots: ArrayLike = None,
    coefs: ArrayLike = None,
    params: dict = {},
) -> dict:
    """Translates a scipy B-spline object or B-spline coefficients and knots for either a iota or pressure profile into a dictionary entries
    that can be handed to `adapt_parameter_file`.

    Args:
        name (str): profile identifyer, has to be either `iota` or `pres`.
        bspl (scipy.interpolate.BSpline): scipy BSpline object. If this is not provided `knots` and `coefs` are expected.
        knots (ArrayLike): Knots for the B-splines. Note that repeated edge knots according to the degree are expected.
        coefs (ArrayLike): Coefficients for the B-splines.
        params (dict, optional): Dictionary of gvec input parameters that will be adapted. Defaults to {}.

    Raises:
        ValueError: If `name` is neither `iota` nor `pres`.
        TypeError: If neither `bspl` nor `knots` and `coefs` is provided.

    Returns:
        dict: Dictionary of gvec input parameters
    """
    if name not in ["iota", "pres"]:
        raise ValueError(
            "Specified profile is not known!"
            + "`which_profile` has to be either `iota` or `pres`."
        )
    if (bspl is None) and (knots is None or coefs is None):
        raise TypeError(
            "`bspl` and at least one of `knots` or `coefs` are None."
            + "Please provide either `bspl` or `knots` and `coefs`"
        )

    if bspl is not None:
        params[f"{name}_coefs"] = bspl.c
        params[f"{name}_knots"] = bspl.t
    else:
        params[f"{name}_coefs"] = coefs
        params[f"{name}_knots"] = knots
    params[f"{name}_type"] = "bspline"

    return params


def logging_setup():
    """Setup default logging configuration for GVEC."""
    import logging

    logging.basicConfig(
        format="{levelname:7s} {message}",
        style="{",
        level=logging.WARNING,
    )
    logging.captureWarnings(True)


def compute_FD(f: np.ndarray, pos, coefs, axis=0):
    """
    1D Finite difference of a function f on equispaced n-dimnesional grid, using FD coefficients coefficients `coefs` and relative integer positions to the central evaluation point `pos`, along one given axis.

    WARNING:
        - if data is periodic, meaning that endpoints of periodic interval are excluded, result can be on all points.
        - If data is not periodic, the result at the **boundaries is WRONG**, for the points |min(pos)| on the left and max(pos) on the right along the given axis.
    Inputs:
        f     : function values on equispaced n-dimnesional grid
        pos   : relative integer positions to the central evaluation point,, as 1d list or 1d array
        coefs : FD coefficients for each position, as 1d list or 1d array, same size as pos!
        axis  : axis along which the FD is computed, default is 0
    Returns:
        df    : Finite-Difference result, same shape as f (see warning above!)
    Examples:
    - examples for first derivative of f:
        - 1st order forward FD: `pos=[1,0]; coefs=[-1,1]/(dx)`
        - 2nd order central FD: `pos=[-1,1]; coefs=[-1,1]/(2*dx)`
        - 4th order central FD: `pos=[-2,-1,1,2]; coefs=[1/12,-2/3,2/3,-1/12]/(dx)`
        - 6th order central FD: `pos=[-3,-2,-1,1,2,3]; coefs=[-1/60,3/20,-3/4,3/4,-3/20,1/60]/(dx)`
        - 8th order central FD: `pos=[-4,-3,-2,-1,1,2,3,4]; coefs=[1/280,-4/105,1/5,-4/5,4/5,-1/5,4/105,-1/280]/(dx)`
    - examples for second derivatives of f:
        - 2nd order central FD: `pos=[-1,0,1]; coefs=[1,-2,1]/(dx**2)`
        - 4th order central FD: `pos=[-2,-1,1,2]; coefs=[-1/12, 4/3,-5/2, 4/3,-1/12]/(dx**2)`
        - 6th order central FD: `pos=[-3,-2,-1,1,2,3]; coefs=[1/90,-3/20,3/2,-49/18,3/2,-3/20,1/90]/(dx**2)`
        - 8th order central FD: `pos=[-4,-3,-2,-1,1,2,3,4]; coefs=[-1/560,8/315,-1/5,8/5,-205/72,8/5,-1/5,8/315,-1/560]/(dx**2)`

    """
    assert axis < f.ndim, f"array does not have the requested dimension {axis}"
    assert len(pos) == len(coefs), (
        f"pos and coefs must have the same length, got {pos} and {coefs}"
    )
    df = np.roll(f, -pos[0], axis=axis) * coefs[0]
    for roll, c in zip(pos[1:], coefs[1:]):
        df += np.roll(f, -roll, axis=axis) * c
    return df


def ellipse_circumference_factor(epsilon: float) -> float:
    """
    Compute the circumference factor of an ellipse with elongation epsilon.
    This uses the approximation by Ramanujan, accurate up to h^5 (6.5% error at ε=10).

    A = a b π = aeff^2 π
    ε = a / b >= 1
    C = 2 π aeff Cf(ε)
    Cf ~ (1 + ε) / (2 √ε) [1 + 3 h / (10 + √(4 - 3 h))]
    h = (ε - 1)^2 / (ε + 1)^2
    """
    h = (epsilon - 1) ** 2 / (epsilon + 1) ** 2
    Cf = (1 + epsilon) / (2 * np.sqrt(epsilon)) * (1 + 3 * h / (10 + np.sqrt(4 - 3 * h)))
    return Cf


def boundary_generator_cases():
    return {
        "ellip_cyl": "elliptic cross-section with no change in zeta",
        "ellip_cyl_breathe": "elliptic cross-section with cross-section area changing with zeta",
        "ellip_cyl_rot": "elliptic cross section of constant ellipticity that only rotates with zeta",
        "ellip_cyl_rot2": "cross section of constant ellipticity that only rotates with zeta, but theta=0 origin remains on positive X1 direction",
        "ellip_cyl_stretch": "elliptic cross section of  changing ellipticity with zeta, not orientation of theta=0",
        "ellip_cyl_helix": "cross section of constant ellipticity, where the axis/boundary moves like a helix over zeta",
        "ellip_cyl_helix_rot": "cross section of constant ellipticity, where the axis/boundary moves like a helix, and rotates along it, over zeta",
        "ellip_cyl_helix_rot2": "cross section of constant ellipticity, where the axis/boundary moves like a helix, and rotates along it, over zeta, but theta=0 origin remains on positive X1 direction",
    }


def boundary_generator(case: str, X1_00=1.0, a0=0.5, ellipticity=0.4, helix_r=0.5):
    """
    Define parameters for some simple boundaries for testing.

    Parameters
    ----------
    case : str
        The name of the boundary:
            see `boundary_generator_cases` dictionary

    X1_00 : float, optional
        =major radius if $X^1=R$
    a0 : float, optional
        =minor radius scale
    ellipticity : float, optional
        =ellipticity of the cross section

    Returns
    -------
    parameter dictionary describing $X^1$ and $X^2$
    """
    params = {}
    match case:
        case "ellip_cyl":
            params["X1_mn_max"] = [1, 0]
            params["X2_mn_max"] = [1, 0]
            params["X1_a_cos"] = {(0, 0): X1_00}
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (1, 0): a0 * (1.0 - ellipticity),
            }
            params["X2_b_sin"] = {
                (1, 0): a0 * (1.0 + ellipticity),
            }
        case "ellip_cyl_breathe":
            breathe = 0.1
            params["X1_mn_max"] = [1, 1]
            params["X2_mn_max"] = [1, 1]
            params["X1_a_cos"] = {(0, 0): X1_00}
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (1, 0): a0 * (1.0 - ellipticity) * (1 + breathe),
                (1, -1): -a0 * (1.0 - ellipticity) * 0.5 * breathe,
                (1, 1): -a0 * (1.0 - ellipticity) * 0.5 * breathe,
            }
            params["X2_b_sin"] = {
                (1, 0): a0 * (1.0 + ellipticity) * (1 + breathe),
                (1, -1): -a0 * (1.0 + ellipticity) * 0.5 * breathe,
                (1, 1): -a0 * (1.0 + ellipticity) * 0.5 * breathe,
            }
        case "ellip_cyl_rot":
            params["X1_mn_max"] = [1, 1]
            params["X2_mn_max"] = [1, 1]
            params["X1_a_cos"] = {(0, 0): X1_00}
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (1, 1): a0,
                (1, -1): -a0 * ellipticity,
            }
            params["X2_b_sin"] = {
                (1, 1): a0,
                (1, -1): a0 * ellipticity,
            }
        case "ellip_cyl_rot2":
            params["X1_mn_max"] = [1, 2]
            params["X2_mn_max"] = [1, 2]
            params["X1_a_cos"] = {(0, 0): X1_00}
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (1, 0): a0,
                (1, -2): -a0 * ellipticity,
            }
            params["X2_b_sin"] = {
                (1, 0): a0,
                (1, -2): a0 * ellipticity,
            }
        case "ellip_cyl_stretch":
            params["X1_mn_max"] = [1, 1]
            params["X2_mn_max"] = [1, 1]
            params["X1_a_cos"] = {(0, 0): X1_00}
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (1, 0): a0,
                (1, 1): -0.5 * a0 * ellipticity,
                (1, -1): -0.5 * a0 * ellipticity,
            }
            params["X2_b_sin"] = {
                (1, 0): a0,
                (1, 1): 0.5 * a0 * ellipticity,
                (1, -1): 0.5 * a0 * ellipticity,
            }

        case "ellip_cyl_helix":
            params["X1_mn_max"] = [1, 1]
            params["X2_mn_max"] = [1, 1]
            params["X1_a_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X2_a_sin"] = {
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
                (1, 0): a0 * (1.0 - ellipticity),
            }
            params["X2_b_sin"] = {
                (1, 0): a0 * (1.0 + ellipticity),
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
        case "ellip_cyl_helix_rot":
            params["X1_mn_max"] = [1, 1]
            params["X2_mn_max"] = [1, 1]
            params["X1_a_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X2_a_sin"] = {
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
                (1, 1): a0,
                (1, -1): -a0 * ellipticity,
            }
            params["X2_b_sin"] = {
                (1, 1): a0,
                (1, -1): a0 * ellipticity,
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
        case "ellip_cyl_helix_rot2":
            params["X1_mn_max"] = [1, 2]
            params["X2_mn_max"] = [1, 2]
            params["X1_a_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X2_a_sin"] = {
                (0, 1): helix_r,
                (0, -1): helix_r,
            }
            params["X1_b_cos"] = {
                (0, 0): X1_00,
                (0, 1): helix_r,
                (0, -1): helix_r,
                (1, 0): a0,
                (1, -2): -a0 * ellipticity,
            }
            params["X2_b_sin"] = {
                (1, 1): a0,
                (1, -1): a0 * ellipticity,
                (0, 0): helix_r,
                (0, -2): helix_r,
            }
        case _:
            raise ValueError(f"request boundary '{case}', does not exist!")

    return params
