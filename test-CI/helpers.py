import os
import numpy as np
from contextlib import contextmanager
from pathlib import Path
import logging
import re
import shutil


# === ASSERTION FUNCTIONS === #


def assert_equal_statefiles(
    pathA: str | Path,
    pathB: str | Path,
    rtol: float = 1e-8,
    atol: float = 1e-12,
):
    """
    Asserts line-wise equality between two files at `pathA` and `pathB`.

    Compare two files line by line and raises an AssertionError if any line in the two files differs.
    Comments (lines starting with '#') are compared exactly, while numerical lines are compared up to the specified
    relative tolerance (`rtol`) and absolute tolerance (`atol`).

    Args:
        pathA (str | Path): The path to the first file.
        pathB (str | Path): The path to the second file.
        rtol (float, optional): The relative tolerance for numerical comparisons. Defaults to 1e-8.
        atol (float, optional): The absolute tolerance for numerical comparisons. Defaults to 1e-12.

    Raises:
        AssertionError: If a line in the two files differs.
    """
    with open(pathA) as fileA, open(pathB) as fileB:
        for lidx, (lineA, lineB) in enumerate(zip(fileA, fileB, strict=True)):
            try:
                # compare comments (exactly)
                if lineA[0] == "#":
                    assert lineA == lineB
                # compare numerically (up to tolerances rtol & atol)
                else:
                    assert np.allclose(
                        *[np.array(l.split(","), float) for l in [lineA, lineB]],
                        rtol=rtol,
                        atol=atol,
                    )
            except AssertionError as e:
                logging.error(f"Line no. {lidx+1} differs for {pathA} and {pathB}")
                raise e


def assert_empty_stderr(path: str | Path = "stderr"):
    """
    Asserts that the specified file (default `stderr`) is empty.

    Args:
        path (str | Path, optional): The path to the stderr file. Defaults to `stderr`.
    """
    with open(path) as file:
        lines = file.readlines()
        assert len(lines) == 0


def assert_stdout_finished(path: str | Path = "stdout"):
    """
    Asserts that the specified file (default `stdout`) ends with "GVEC SUCESSFULLY FINISHED!".

    Args:
        path (str | Path, optional): The path to the stdout file. Defaults to `stdout`.
    """
    with open(path) as file:
        lines = file.readlines()
        assert lines[-2].startswith(" GVEC SUCESSFULLY FINISHED!")


# === HELPER FUNCTIONS === #


@contextmanager
def chdir(target: str | Path):
    """
    Contextmanager to change working directory, returns to the original directory when exiting.

    Args:
        target (str | Path): The target directory to change to.

    Example:
    >>> with chdir('/path/to/directory'):
    ...     # Code executed within the new directory
    ... # Code executed in the original directory
    """
    source = os.getcwd()
    try:
        os.chdir(target)
        yield
    finally:
        os.chdir(source)


def adapt_parameter_file(source: str | Path, target: str | Path, **kwargs):
    """
    Copy the `source` file to the `target` file, replacing the parameters according to `kwargs`.

    Args:
        source (str or Path): The path to the source parameter file.
        target (str or Path): The path to the target parameter file.
        **kwargs: Keyword arguments representing the parameters to be replaced.

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
    with open(source, "r") as source_file, open(target, "w") as target_file:
        for line in source_file:
            if m := re.match(
                r"([^!]*)("
                + "|".join(kwargs.keys())
                + r")(\s*=\s*)(\(.*\)|[^!\s]*)(.*)",
                line,
            ):
                prefix, key, sep, value, suffix = m.groups()
                line = f"{prefix}{key}{sep}{kwargs[key]}{suffix}\n"
                occurrences[key] += 1
            target_file.write(line)
        # add key,value pair if not existing in parameterfile.
        for key, v in occurrences.items():
            if v == 0:
                target_file.write(f"\n {key} = {kwargs[key]}")
                occurrences[key] += 1
    assert all(
        [v == 1 for v in occurrences.values()]
    ), f"bad number of occurrences: {occurrences}"
