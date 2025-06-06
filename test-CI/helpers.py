import os
import numpy as np
from contextlib import contextmanager
from pathlib import Path
import logging
import re
import shutil
import sys

# === ASSERTION FUNCTIONS === #


def check_diff_files(
    pathA: str | Path,
    pathB: str | Path,
    rtol: float = 1e-8,
    atol: float = 1e-12,
    ignore_lines: list[int] = [],
    ignore_columns: list[int] = [],
    ignore_regexs: list[str] = [],
    warn_regexs: list[str] = [],
    logger: logging.Logger = None,
) -> tuple[int, int]:
    """
    Asserts line-wise equality between two files at `pathA` and `pathB`.

    Compare two files line by line and raises an AssertionError if any line in the two files differs.
    Float values (of the form 1.234E+45 or -1.234E-45) are compared up to the specified relative and absolute tolerance (`rtol`) and (`atol`),
    while all other content is compared exactly (leading and trailing whitespace is ignored).

    Args:
        pathA (str | Path): The path to the first file.
        pathB (str | Path): The path to the second file.
        rtol (float, optional): The relative tolerance for numerical comparisons. Defaults to 1e-8.
        atol (float, optional): The absolute tolerance for numerical comparisons. Defaults to 1e-12.
        ignore_lines (list[int], optional): List of line indices to ignore during comparison. Defaults to [].
        ignore_columns (list[int], optional): List of column indices to ignore during numerical comparison. Defaults to [].
        ignore_regexs (list[str], optional): List of regular expressions to match lines to ignore during comparison. Defaults to [].
        warn_regexs (list[str], optional): List of regular expressions to match lines to suppress and warn about during comparison. Defaults to [].
        logger (logging.Logger, optional): The logger to use for logging. If not provided, falls back to the root logger.

    Raises:
        AssertionError: If a line in the two files differs.

    Returns:
        tuple[int, int]: The number of textual differences and the number of numerical differences.
    """
    # fall back to root logger
    if logger is None:
        logger = logging.getLogger()
    # count the number of differences
    txt_differences = 0
    num_differences = 0
    warnings = 0
    # read files & filter lines based on ignore_lines and ignore_regexs
    linesA, linesB, lidxsA, lidxsB = [], [], [], []
    for lines, lidxs, path in [(linesA, lidxsA, pathA), (linesB, lidxsB, pathB)]:
        with open(path) as file:
            for lidx, line in enumerate(file):
                if lidx in ignore_lines or any(
                    [re.search(regex, line) for regex in ignore_regexs]
                ):
                    continue
                lines.append(line)
                lidxs.append(lidx)

    # filter lines based on warn_regexs
    lineZA, lineZB = zip(lidxsA, linesA), zip(lidxsB, linesB)
    warningsA, warningsB = [], []
    for (lidxA, lineA), (lidxB, lineB) in zip(lineZA, lineZB):
        while True:
            if any(
                [
                    re.search(regex, lineA) and not re.search(regex, lineB)
                    for regex in warn_regexs
                ]
            ):
                warnings += 1
                logger.warning(f"Extra lineA #{lidxA + 1} (suppressed)")
                logger.debug(f"=> Line A: {lineA!r}")
                warningsA.append(lidxA)
                try:
                    lidxA, lineA = next(lineZA)
                except StopIteration:
                    lidxA, lineA = None, ""
            elif any(
                [
                    re.search(regex, lineB) and not re.search(regex, lineA)
                    for regex in warn_regexs
                ]
            ):
                warnings += 1
                logger.warning(f"Extra lineB #{lidxB + 1} (suppressed)")
                logger.debug(f"=> Line B: {lineB!r}")
                warningsB.append(lidxB)
                try:
                    lidxB, lineB = next(lineZB)
                except StopIteration:
                    lidxB, lineB = None, ""
            else:
                break
    if len(warningsA) > 0:
        lidxsA, linesA = zip(
            *[
                (lidx, line)
                for lidx, line in zip(lidxsA, linesA)
                if lidx not in warningsA
            ]
        )
    if len(warningsB) > 0:
        lidxsB, linesB = zip(
            *[
                (lidx, line)
                for lidx, line in zip(lidxsB, linesB)
                if lidx not in warningsB
            ]
        )
    # compare number of lines
    txt_differences += abs(len(linesA) - len(linesB))
    if txt_differences > 0:
        if Path(pathA).name != Path(pathB).name:
            logger.info(f"--- Comparing {Path(pathA).name} and {Path(pathB).name} ---")
        else:
            logger.info(f"--- Comparing {Path(pathA).name} ---")
        logger.error(
            f"files (after applying filters) do not have the same number of lines ({len(linesA)} vs {len(linesB)})"
        )
        return txt_differences, num_differences, warnings
    # compare lines
    for lidxA, lidxB, lineA, lineB in zip(lidxsA, lidxsB, linesA, linesB, strict=True):
        # compare with regex
        # split line into text and float parts (where floats have the format 1.234E+45 or -1.234E-45)
        splitA, splitB = (
            re.split(r"(-?\d\.\d+E[+-]\d+)", line) for line in (lineA, lineB)
        )
        # extract the text fragments
        textsA, textsB = (
            [text.strip() for text in split[::2]] for split in (splitA, splitB)
        )
        if textsA != textsB:
            if num_differences == 0:
                if Path(pathA).name != Path(pathB).name:
                    logger.info(
                        f"--- Comparing {Path(pathA).name} and {Path(pathB).name} ---"
                    )
                else:
                    logger.info(f"--- Comparing {Path(pathA).name} ---")
            if not any([re.search(regex, lineA) for regex in warn_regexs]):
                txt_differences += 1
                logger.error(
                    f"LineA #{lidxA + 1} and LineB #{lidxB + 1} differ textually"
                )
            else:
                warnings += 1
                logger.warning(
                    f"LineA #{lidxA + 1} and LineB #{lidxB + 1} differ textually (suppressed)"
                )
            logger.debug(f"=> Line A: {lineA!r}")
            logger.debug(f"=> Line B: {lineB!r}")
        # extract the float fragments
        floatsA, floatsB = (
            np.array(split[1::2], dtype=float) for split in (splitA, splitB)
        )
        if not (len(floatsA) == len(floatsB)):
            num_differences += 1
            logger.error(
                f"LineA #{lidxA + 1} and LineB #{lidxB + 1} do not have the same number of floats"
            )
        else:
            select = np.ones(floatsA.shape, dtype=bool)
            select[[i for i in ignore_columns if i < select.size]] = False
            floatsA, floatsB = floatsA[select], floatsB[select]
            close = np.isclose(floatsA, floatsB, rtol=rtol, atol=atol)
            # error pattern for difference
            pattern = "".join("." if c else "x" for c in close)
            # assert equality with better error messages
            if not all(close):
                if num_differences == 0:
                    if Path(pathA).name != Path(pathB).name:
                        logger.info(
                            f"--- Comparing {Path(pathA).name} and {Path(pathB).name} ---"
                        )
                    else:
                        logger.info(f"--- Comparing {Path(pathA).name} ---")
                if not (any([re.search(regex, lineA) for regex in warn_regexs])):
                    num_differences += 1
                    logger.error(
                        f"LineA #{lidxA + 1} and LineB #{lidxB + 1} differ  numerically (rtol={rtol}, atol={atol}, {pattern})"
                    )
                else:
                    logger.warning(
                        f"LineA #{lidxA + 1} and LineB #{lidxB + 1} differ numerically (rtol={rtol}, atol={atol}, {pattern}) (suppressed)"
                    )
                logger.debug(f"=> Line A: {lineA!r}")
                logger.debug(f"=> Line B: {lineB!r}")
                with np.printoptions(formatter={"float": "{:.2E}".format}):
                    logger.debug(f"=> |A-B|: {abs(floatsA - floatsB)}")
    return txt_differences, num_differences, warnings


def assert_empty_stderr(path: str | Path = "stderr.txt", slurm: bool = False):
    """
    Asserts that the specified file (default `stderr.txt`) is empty.

    Args:
        path (str | Path, optional): The path to the stderr file. Defaults to `stderr.txt`.

    """
    with open(path) as file:
        lines = file.readlines()
    if slurm is True:
        # check for SLURM stderr in stderr.txt
        assert len(lines) == 2, "Errors found in stderr.txt"
        assert ("srun:" in line for line in lines), (
            "SLURM stderr not found in stderr.txt"
        )
    else:
        # check for an empty stderr.txt file
        assert len(lines) == 0


def assert_stdout_finished(
    path: str | Path = "stdout.txt", message="SUCCESSFULLY FINISHED"
):
    """
    Asserts that the specified file (default `stdout.txt`) ends with a message in the second to last line.

    Args:
        path (str | Path, optional): The path to the stdout.txt file. Defaults to `stdout.txt`.
    """
    with open(path) as file:
        lines = file.readlines()
        assert message.strip() in lines[-2], (
            f"final message '{message.strip()}' not found in stdout.txt"
        )


def assert_stdout_OpenMP_MPI(path: str | Path = "stdout.txt"):
    with open(path) as file:
        lines = file.readlines()
    # check OpenMP
    if os.environ.get("OMP_MODE") == "ompON":
        assert any(re.search(r"OpenMP threads", line) for line in lines), (
            "no OpenMP support despite ompON"
        )
    elif os.environ.get("OMP_MODE") == "ompOFF":
        assert not any(re.search(r"OpenMP threads", line) for line in lines), (
            "OpenMP support despite ompOFF"
        )
    # check MPI
    if os.environ.get("MPI_MODE") == "mpiON":
        assert any(re.search(r"MPI tasks", line) for line in lines), (
            "no MPI support despite mpiON"
        )
    elif os.environ.get("MPI_MODE") == "mpiOFF":
        assert not any(re.search(r"MPI tasks", line) for line in lines), (
            "MPI support despite mpiOFF"
        )


def assert_stdout_no_NaN(path: str | Path = "stdout.txt"):
    with open(path) as file:
        lines = file.readlines()
        assert not any(re.search(r"NaN", line) for line in lines), "found NaN in stdout"


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
