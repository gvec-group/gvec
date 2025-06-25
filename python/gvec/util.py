# Copyright (c) 2025 GVEC Contributors, Max Planck Institute for Plasma Physics
# License: MIT
"""GVEC utility module

This module is part of the gvec python package, but also used directly in the tests.
"""

import contextlib
import os
import re
import shutil
from collections.abc import Mapping, MutableMapping, Iterable
from pathlib import Path
from typing import Literal
from copy import deepcopy

from numpy.typing import ArrayLike

try:
    from scipy.interpolate import BSpline
except ImportError:
    BSpline = None


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
    FLOAT = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
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


def flip_parameters_theta(parameters: MutableMapping) -> MutableMapping:
    import copy

    parameters2 = copy.deepcopy(parameters)
    if "X1_b_cos" in parameters:
        for (m, n), value in parameters["X1_b_cos"].items():
            if m == 0:
                continue
            parameters2["X1_b_cos"][m, -n] = value
    if "X1_b_sin" in parameters:
        for (m, n), value in parameters["X1_b_sin"].items():
            if m == 0:
                continue
            parameters2["X1_b_sin"][m, -n] = -value
    if "X2_b_cos" in parameters:
        for (m, n), value in parameters["X2_b_cos"].items():
            if m == 0:
                continue
            parameters2["X2_b_cos"][m, -n] = value
    if "X2_b_sin" in parameters:
        for (m, n), value in parameters["X2_b_sin"].items():
            if m == 0:
                continue
            parameters2["X2_b_sin"][m, -n] = -value
    return parameters2


def flip_parameters_zeta(parameters: MutableMapping) -> MutableMapping:
    import copy

    parameters2 = copy.deepcopy(parameters)
    for var in ["X1_b", "X2_b"]:
        if f"{var}_cos" in parameters:
            for (m, n), value in parameters[f"{var}_cos"].items():
                if m == 0:
                    continue
                parameters2[f"{var}_cos"][m, -n] = value
        if f"{var}_sin" in parameters:
            for (m, n), value in parameters[f"{var}_sin"].items():
                if m == 0:
                    parameters2[f"{var}_sin"][m, n] = -value
                else:
                    parameters2[f"{var}_sin"][m, -n] = value
    for var in ["X1_a", "X2_a"]:
        if f"{var}_sin" in parameters:
            for (m, n), value in parameters[f"{var}_sin"].items():
                assert m == 0
                parameters2[f"{var}_sin"][m, n] = -value
        if f"{var}_cos" in parameters:
            for (m, n), value in parameters[f"{var}_cos"].items():
                assert m == 0
                # parameters2[f"{var}_cos"][m, n] = value
    return parameters2


def parameters_from_vmec(nml: Mapping) -> CaseInsensitiveDict:
    import numpy as np

    M, N = nml["mpol"] - 1, nml["ntor"]
    stellsym = nml["lasym"]  # stellarator symmetry
    params = CaseInsensitiveDict(
        {
            "nfp": nml["nfp"],
            "X1_mn_max": f"(/{M}, {N}/)",
            "X2_mn_max": f"(/{M}, {N}/)",
            "LA_mn_max": f"(/{M}, {N}/)",
            "PHIEDGE": nml["phiedge"],
        }
    )
    if stellsym:
        params["X1_sin_cos"] = "_cos_"
        params["X2_sin_cos"] = "_sin_"
        params["LA_sin_cos"] = "_sin_"
    else:
        params["X1_sin_cos"] = "_sincos_"
        params["X2_sin_cos"] = "_sincos_"
        params["LA_sin_cos"] = "_sincos_"

    # --- boundary --- #
    rbc = np.array(nml["rbc"], dtype=float)
    zbs = np.array(nml["zbs"], dtype=float)
    if not rbc.shape == zbs.shape == (M + 1, 2 * N + 1):
        raise ValueError(
            f"VMEC namelist arrays 'rbc' and 'zbs' have shape {rbc.shape} and {zbs.shape} that does not match the expected shape {(M + 1, 2 * N + 1)=}"
        )
    if not stellsym:
        rbs = np.array(nml["rbs"], dtype=float)
        zbc = np.array(nml["zbc"], dtype=float)
        if not rbs.shape == zbc.shape == (M + 1, 2 * N + 1):
            raise ValueError(
                f"VMEC namelist arrays 'rbs' and 'zbc' have shape {rbs.shape} and {zbc.shape} that does not match the expected shape {(M + 1, 2 * N + 1)=}"
            )

    params["X1_b_cos"] = {}
    params["X2_b_sin"] = {}
    if not stellsym:
        params["X1_b_sin"] = {}
        params["X2_b_cos"] = {}
    for m in range(M + 1):
        for n in range(-N, N + 1):
            if m == 0 and n < 0:
                continue
            params["X1_b_cos"][m, n] = rbc[m, n + N]
            if not stellsym:
                params["X1_b_sin"][m, n] = rbs[m, n + N]
                params["X2_b_cos"][m, n] = zbc[m, n + N]
            params["X2_b_sin"][m, n] = zbs[m, n + N]

    # --- axis --- #
    params["X1_a_cos"] = {(0, n): v for n, v in enumerate(nml["raxis_cc"])}
    params["X2_a_sin"] = {(0, n): v for n, v in enumerate(nml["zaxis_cs"])}
    if not stellsym and nml["raxis_cs"] is not None:
        params["X1_a_sin"] = {(0, n): v for n, v in enumerate(nml["raxis_cs"])}
    if not stellsym and nml["zaxis_cc"] is not None:
        params["X2_a_cos"] = {(0, n): v for n, v in enumerate(nml["zaxis_cc"])}

    return params


def axis_from_boundary(parameters: MutableMapping) -> MutableMapping:
    import copy

    parameters2 = copy.deepcopy(parameters)
    N = parameters["X1_mn_max"][1]
    parameters2["X1_a_cos"] = {parameters["X1_b_cos"][0, n] for n in range(N + 1)}
    parameters2["X2_a_sin"] = {parameters["X2_b_sin"][0, n] for n in range(N + 1)}
    if "X1_b_sin" in parameters:
        parameters2["X1_a_sin"] = {parameters["X1_b_sin"][0, n] for n in range(N + 1)}
    if "X2_b_cos" in parameters:
        parameters2["X2_a_cos"] = {parameters["X2_b_cos"][0, n] for n in range(N + 1)}
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
                output[key][f"({m}, {n:2d})"] = val
        elif key.lower() == "stages":
            output[key] = [stringify_mn_parameters(stage) for stage in value]
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
