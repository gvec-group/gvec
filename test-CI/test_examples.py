import pytest
import sys, os
import subprocess
import re
from pathlib import Path
import shutil

import helpers


# === HELPER FUNCTIONS === #


def discover_subdirs(path: Path | str = ".") -> list[str]:
    """
    Discover subdirectories in this module's directory for dynamic test detection.

    Ignores directories starting with `_` or `.`.

    Args:
        path (optional): Path to the directory to discover subdirectories in.\
            Defaults to the directory of this file.\
            Relative paths are relative to the directory of this file.

    Returns:
        list: A sorted list of subdirectory names.
    """
    # append relative path to the path of this file
    path = Path(path)
    if not path.is_absolute():
        path = Path(__file__).parent / path
    assert path.is_dir()
    # list subdirectories except for those starting with `_` or `.`
    return sorted(
        [
            sd
            for sd in os.listdir(path)
            if (path / sd).is_dir()
            and not sd.startswith("_")
            and not sd.startswith(".")
        ]
    )


def discover_examples(generate_shortruns: bool = True):
    # discover example testcases
    examples = discover_subdirs("examples")
    # add example marker
    testcases = [pytest.param("examples", example, marks=pytest.mark.example) for example in examples]
    if not generate_shortruns:
        return testcases
    # generate shortruns
    with helpers.chdir(Path(__file__).parent):
        directory = Path("shortruns")
        if directory.exists():
            shutil.rmtree(directory)
        # create a directory for each example, copy `parameters.ini` and modify the `testlevel`
        directory.mkdir()
        for example in examples:
            subdir = directory / example
            subdir.mkdir()
            helpers.adapt_parameter_file("examples" / Path(example) / "parameter.ini", subdir / "parameter.ini", testlevel=2)
            testcases.append(pytest.param("shortruns", example, marks=pytest.mark.shortrun))
    return testcases


# === TESTS === #


@pytest.mark.parametrize(["stage", "testcase"], discover_examples())
def test_examples(binpath, stage, testcase):
    """test end2end GVEC run with `{stage}/{testcase}/parameter.ini`"""
    directory = Path(__file__).parent / stage / testcase
    args = [binpath, "parameter.ini"]
    # find restart statefile
    if "restart" in testcase:
        base_directory = Path(__file__).parent / stage / testcase[:-8]  # [-8] to remove `_restart` suffix
        states = [sd for sd in os.listdir(base_directory) if "State" in sd and sd.endswith(".dat")]
        args.append(base_directory / sorted(states)[-1])
    # run gvec
    with helpers.chdir(directory):
        # run GVEC
        with open("stdout", "w") as stdout, open("stderr", "w") as stderr:
            subprocess.run(args, text=True, stdout=stdout, stderr=stderr)
        # check if GVEC was successful
        helpers.assert_empty_stderr()
        helpers.assert_stdout_finished()


@pytest.mark.regression
@pytest.mark.parametrize("testcase", discover_subdirs("examples"))
def test_regression(rootdir, refdir, testcase, logger):
    """test regression of example GVEC runs and restarts"""
    directory = Path(__file__).parent / "examples" / testcase
    refdirectory = refdir / directory.relative_to(rootdir)
    assert refdirectory.exists(), f"Reference testcase {refdirectory} does not exist"
    # compare output to reference
    for filename in os.listdir(directory):
        if re.match(r'\w+State\w+\.dat', filename) and filename in os.listdir(refdirectory):
            helpers.assert_equal_statefiles(
                directory / filename,
                refdirectory / filename,
            )
