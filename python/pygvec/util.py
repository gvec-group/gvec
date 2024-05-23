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

from pathlib import Path
import re
import shutil


def adapt_parameter_file(source: str | Path, target: str | Path, **kwargs):
    """
    Copy the `source` file to the `target` file and replace the parameters according to `kwargs`.

    Args:
        source (str or Path): The path to the source parameter file.
        target (str or Path): The path to the target parameter file.
        **kwargs: Keyword arguments representing the parameters to be replaced.
                  if the value of the key is "!", the line with the keyword is uncommented (must exist in the parameterfile)

    Raises:
        AssertionError: If the number of occurrences for any parameter is not exactly 1.

    Notes:
        - If no parameters are provided in `kwargs`, the function simply copies the `source` file to the `target` file.
        - The function replaces the parameters in the format `key = value`, where value is either a sequence of characters containing
        no whitespace or a single pair of parentheses with any content. The value from `kwargs` is inserted using the standard python
        string conversion. There may be a comment, starting with `!`, after the value.
        - If a parameter already exists in the `source` file, its value is replaced with the corresponding value from `kwargs`.
        - If a parameter does not exist in the `source` file, it is added to the `target` file.

    Example:
    >>> adapt_parameter_file('/path/to/source.ini', '/path/to/target.ini', param1=1.2, param2="(1, 2, 3)")
    """
    if not len(kwargs.keys()):
        shutil.copy2(source, target)
        return
    occurrences = {key: 0 for key in kwargs}
    assert target != source
    with open(source, "r") as source_file, open(target, "w") as target_file:
        for line in source_file:
            if m := re.match(
                r"\s*([^!=\s\(]+)\s*\(\s*([-\d]+);\s*([-\d]+)\)\s*=\s*([-+\d\.Ee]+)",
                line,
            ):
                key, *mn, value = m.groups()
                if key in kwargs:
                    if (int(mn[0]), int(mn[1])) in kwargs[key]:
                        line = f"{key}({mn[0]};{mn[1]}) = {kwargs[key][(int(mn[0]), int(mn[1]))]}\n"
                        occurrences[key, int(mn[0]), int(mn[1])] += 1
            elif m := re.match(
                r"([\s!]*)("
                + "|".join(
                    [
                        key
                        for key, value in kwargs.items()
                        if not isinstance(value, dict)
                    ]
                )
                + r")(\s*=\s*)(\([^\)]*\)|[^!\s]*)(.*)",
                line,
            ):
                prefix, key, sep, value, suffix = m.groups()
                if "!" in prefix:  # found commented keyword
                    if kwargs[key] == "!":  # only uncomment keyword
                        line = f"{key}{sep}{value}{suffix}\n"
                        occurrences[key] += 1
                else:  # found uncommented keywords
                    if not (kwargs[key] == "!"):  # use new keyword
                        line = f"{prefix}{key}{sep}{kwargs[key]}{suffix}\n"
                        occurrences[key] += 1
                    else:  # use the existing keyword,value pair with a comment
                        line = f"{prefix}{key}{sep}{value} !!WAS ALREADY UNCOMMENTED!! {suffix}\n"
                        occurrences[key] += 1
            target_file.write(line)
        # add key,value pair if not existing in parameterfile.
        for key, v in occurrences.items():
            if v == 0:
                if not (
                    kwargs[key] == "!"
                ):  # ignore uncommenting keywords if not found
                    if isinstance(kwargs[key], dict):
                        for (m, n), value in kwargs[key].items():
                            target_file.write(f"\n {key}({m};{n}) = {value}")
                    else:
                        target_file.write(f"\n {key} = {kwargs[key]}")
                    occurrences[key] += 1
    assert all(
        [v == 1 for v in occurrences.values()]
    ), f"bad number of occurrences in adapt_parameter_file: {occurrences}"


def read_parameter_file(path: str | Path) -> dict:
    """
    Read the parameters from the specified parameter file.

    Args:
        path (str | Path): The path to the parameter file.

    Returns:
        dict: A dictionary containing the parameters from the parameter file.

    Example:
    >>> read_parameter_file('/path/to/parameter.ini')
    {'param1': 1.2, 'param2': (1, 2, 3), 'param3': {(-1, 0): 0.5, (0, 0): 1.0}}
    """
    parameters = {}
    with open(path, "r") as file:
        for line in file:
            # match parameter in the form `key(m,n) = value` with m,n integers and value a float
            if ln := re.match(
                r"\s*([^!=\s\(]+)\s*\(\s*([-\d]+);\s*([-\d]+)\)\s*=\s*([-+\d\.Ee]+)",
                line,
            ):
                key, m, n, value = ln.groups()
                if key not in parameters:
                    parameters[key] = {(int(m), int(n)): float(value)}
                else:
                    assert (int(m), int(n)) not in parameters[key]
                    parameters[key][int(m), int(n)] = float(value)
            # match parameter in the form `key = (/value/)` with value an array of ints or floats, separated by commas
            elif ln := re.match(r"\s*([^!=\s]+)\s*=\s*\(/(.*)/\)", line):
                key, value = ln.groups()
                assert key not in parameters
                try:
                    parameters[key] = tuple(map(int, value.split(",")))
                except ValueError:
                    parameters[key] = tuple(map(float, value.split(",")))
            # match parameter in the form `key = value`
            elif ln := re.match(r"\s*([^!=\s]+)\s*=\s*([^!\s]*)", line):
                key, value = ln.groups()
                assert key not in parameters
                try:
                    parameters[key] = int(value)
                except ValueError:
                    try:
                        parameters[key] = float(value)
                    except ValueError:
                        parameters[key] = value
    return parameters
