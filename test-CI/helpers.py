import os
import numpy as np
from contextlib import contextmanager
from pathlib import Path
import logging
import re
import shutil


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
    # compare files
    with open(pathA) as fileA, open(pathB) as fileB:
        # filter lines based on ignore_lines and ignore_regexs
        linesA, linesB, lidxsA, lidxsB = [], [], [], [] 
        for lidxs, lines, file in zip((lidxsA, lidxsB), (linesA, linesB), (fileA, fileB)):
            for lidx, line in enumerate(file):
                if not (lidx in ignore_lines or any([re.match(regex, line) for regex in ignore_regexs])):
                    lines.append(line)
                    lidxs.append(lidx)
        txt_differences += abs(len(linesA) - len(linesB))
        if(txt_differences > 0):
            if Path(pathA).name != Path(pathB).name:
                logger.info(f"--- Comparing {Path(pathA).name} and {Path(pathB).name} ---")
            else:
                logger.info(f"--- Comparing {Path(pathA).name} ---")
            logger.error(f"files (after applying filters) do not have the same number of lines ({len(linesA)} vs {len(linesB)})")
            return txt_differences, num_differences
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
                        logger.info(f"--- Comparing {Path(pathA).name} and {Path(pathB).name} ---")
                    else:
                        logger.info(f"--- Comparing {Path(pathA).name} ---")
                if not (any([re.match(regex, lineA) for regex in warn_regexs])):
                    txt_differences += 1
                    logger.error(f"LineA #{lidxA+1} and LineB #{lidxB+1} differ textually")
                else:
                    logger.warning(f"LineA #{lidxA+1} and LineB #{lidxB+1} differ textually (suppressed)")
                logger.debug(f"=> Line A: {lineA!r}")
                logger.debug(f"=> Line B: {lineB!r}")
            # extract the float fragments
            floatsA, floatsB = (
                np.array(split[1::2], dtype=float) for split in (splitA, splitB)
            )
            if not (len(floatsA)==len(floatsB)):
                num_differences +=1
                logger.error(f"LineA #{lidxA+1} and LineB #{lidxB+1} do not have the same number of floats")
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
                            logger.info(f"--- Comparing {Path(pathA).name} and {Path(pathB).name} ---")
                        else:
                            logger.info(f"--- Comparing {Path(pathA).name} ---")
                    if not (any([re.match(regex, lineA) for regex in warn_regexs])):
                        num_differences += 1
                        logger.error(f"LineA #{lidxA+1} and LineB #{lidxB+1} differ  numerically (rtol={rtol}, atol={atol}, {pattern})")
                    else:
                        logger.warning(f"LineA #{lidxA+1} and LineB #{lidxB+1} differ numerically (rtol={rtol}, atol={atol}, {pattern}) (suppressed)")
                    logger.debug(f"=> Line A: {lineA!r}")
                    logger.debug(f"=> Line B: {lineB!r}")
                    logger.debug(f"=> |A-B|: {abs(floatsA - floatsB)}")
    return txt_differences, num_differences


def assert_empty_stderr(path: str | Path = "stderr.txt"):
    """
    Asserts that the specified file (default `stderr.txt`) is empty.

    Args:
        path (str | Path, optional): The path to the stderr file. Defaults to `stderr.txt`.
    """
    with open(path) as file:
        lines = file.readlines()
        assert len(lines) == 0


def assert_stdout_finished(
    path: str | Path = "stdout.txt", message="SUCESSFULLY FINISHED"
):
    """
    Asserts that the specified file (default `stdout.txt`) ends with a message in the second to last line.

    Args:
        path (str | Path, optional): The path to the stdout.txt file. Defaults to `stdout.txt`.
    """
    with open(path) as file:
        lines = file.readlines()
        assert (
            message.strip() in lines[-2]
        ), f"final message '{message.strip()}' not found in stdout.txt"


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
    with open(source, "r") as source_file, open(target, "w") as target_file:
        for line in source_file:
            if m := re.match(
                r"([\s!]*)("
                + "|".join(kwargs.keys())
                + r")(\s*=\s*)(\(.*\)|[^!\s]*)(.*)",
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
                    target_file.write(f"\n {key} = {kwargs[key]}")
                    occurrences[key] += 1

    assert all(
        [v == 1 for v in occurrences.values()]
    ), f"bad number of occurrences in adapt_parameter_file: {occurrences}"
