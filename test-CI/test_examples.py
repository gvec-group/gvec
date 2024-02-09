import pytest
import sys, os
import subprocess
import re
from pathlib import Path

import helpers

# === AUTOMATIC / SETUP FIXTURES === #


@pytest.fixture(autouse=True, scope="module")
def chdir():
    """change working directory to this module's directory"""
    with helpers.chdir(os.path.dirname(__file__)):
        yield  # ToDo: does this work as expected?


# === HELPER FUNCTIONS === #


def discover_subdirs(path: Path | str = ".", restart: bool | None = None) -> list[str]:
    """
    Discover subdirectories in this module's directory for dynamic test detection.

    Ignores directories starting with `_` or `.`.

    Args:
        path (optional): Path to the directory to discover subdirectories in.\
            Defaults to the directory of this file.\
            Relative paths are relative to the directory of this file.
        restart (optional): Specifies whether to restart the discovery process.\
            None (default) means all directories are discovered.\
            True means only directories with a `_restart` suffix are discovered.\
            False means only directories without a `_restart` suffix are discovered.

    Returns:
        list: A sorted list of subdirectory names.
    """
    # append relative path to the path of this file
    path = Path(path)
    if not path.is_absolute():
        path = Path(__file__).parent / path
    assert path.is_dir()
    # list subdirectories except for those starting with `_` or `.`
    subdirs = sorted(
        [
            sd
            for sd in os.listdir(path)
            if (path / sd).is_dir()
            and not sd.startswith("_")
            and not sd.startswith(".")
        ]
    )
    # filter by `_restart` suffix
    if restart is True:
        return [sd[:-8] for sd in subdirs if sd.endswith("_restart")]
    elif restart is False:
        return [sd for sd in subdirs if not sd.endswith("_restart")]
    else:
        return subdirs


# === TESTS === #


@pytest.mark.example
@pytest.mark.parametrize("testcase", discover_subdirs("examples", restart=False))
def test_examples(binpath, testcase):
    """test end2end GVEC run with `examples/testcase/parameter.ini`"""
    directory = Path(__file__).parent / "examples" / testcase
    with helpers.chdir(directory):
        # run GVEC
        with open("stdout", "w") as stdout, open("stderr", "w") as stderr:
            subprocess.run(
                [binpath, "parameter.ini"], text=True, stdout=stdout, stderr=stderr
            )
        # check if GVEC was successful
        helpers.assert_empty_stderr()
        helpers.assert_stdout_finished()


@pytest.mark.example
@pytest.mark.parametrize("testcase", discover_subdirs("examples", restart=True))
def test_restart(binpath, testcase):
    """test end2end GVEC restart-run with `parameter.ini`"""
    # look for last state in base directory
    base_directory = Path(__file__).parent / "examples" / testcase
    states = [sd for sd in os.listdir(base_directory) if "State" in sd and sd.endswith(".dat")]
    last_state = base_directory / sorted(states)[-1]
    # run GVEC in restart mode
    directory = Path(__file__).parent / "examples" / f"{testcase}_restart"
    with helpers.chdir(directory):
        with open("stdout", "w") as stdout, open("stderr", "w") as stderr:
            subprocess.run(
                [binpath, "parameter.ini", last_state], text=True, stdout=stdout, stderr=stderr
            )
        # check if GVEC was successful
        helpers.assert_empty_stderr()
        helpers.assert_stdout_finished()


@pytest.mark.parametrize("directory", discover_subdirs())
def test_regression_gvec(binpath, rootpath, refpath, directory):
    """test end2end GVEC run (and restarts)"""
    refdirectory = refpath / Path(__file__).parent.relative_to(rootpath) / directory
    assert refdirectory.exists(), f"Reference testcase {refdirectory} does not exist"
    with helpers.chdir(directory):
        # compare output to reference
        for path in os.listdir():
            if match := re.match(r"(REF_(\w+\.dat))", path):
                helpers.assert_equal_statefiles(
                    match.groups()[0],
                    match.groups()[1],
                )
