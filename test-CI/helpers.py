import os
import numpy as np
from contextlib import contextmanager
from pathlib import Path
import logging
import re
import shutil


def assert_equal_statefiles(
    pathA: str | Path,
    pathB: str | Path,
    rtol: float = 1e-8,
    atol: float = 1e-12,
):
    """assert line-wise equality between to files at `pathA` and `pathB`"""
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


def assert_empty_stderr(path="stderr"):
    with open(path) as file:
        lines = file.readlines()
        assert len(lines) == 0


def assert_stdout_finished(path="stdout"):
    with open(path) as file:
        lines = file.readlines()
        assert lines[-2].startswith(" GVEC SUCESSFULLY FINISHED!")


@contextmanager
def chdir(target):
    """context to change working directory, returns to original directory when exiting"""
    source = os.getcwd()
    try:
        os.chdir(target)
        yield
    finally:
        os.chdir(source)


def adapt_parameter_file(source: str | Path, target: str | Path, **kwargs):
    """copy `source` to `target`, replacing the parameters according to `kwargs`"""
    if not len(kwargs.keys()):
        shutil.copy2(source, target)
        return
    occurances = {key: 0 for key in kwargs}
    with open(source, "r") as source_file, open(target, "w") as target_file:
        for line in source_file:
            if m := re.match(r'([^!]*)(' + "|".join(kwargs.keys()) + r')(\s*=\s*)(\(.*\)|[^!\s]*)(.*)', line):
                prefix, key, sep, value, suffix = m.groups()
                line = f"{prefix}{key}{sep}{kwargs[key]}{suffix}\n"
                occurances[key] += 1
            target_file.write(line)
        # add key,value pair if not existing in parameterfile.
        for key, v in occurances.items():
            if v == 0:
                target_file.write(f"\n {key} = {kwargs[key]}")
                occurances[key] += 1       
    assert all([v == 1 for v in occurances.values()]), f"bad number of occurances: {occurances}"
