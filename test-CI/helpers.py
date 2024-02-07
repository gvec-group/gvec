import os
import numpy as np
from contextlib import contextmanager
from typing import Union
import logging


def assert_equal_statefiles(
    pathA: Union[str, bytes, os.PathLike],
    pathB: Union[str, bytes, os.PathLike],
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
