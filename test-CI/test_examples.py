import pytest
import sys, os
import subprocess
import re
from pathlib import Path
import shutil

import helpers


# === HELPER FUNCTIONS & FIXTURES === #


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


@pytest.fixture(scope="session", params=discover_subdirs("examples"))
def testcase(request):
    """
    Fixture to parametrize the `testcase`

    Discoveres subdirectories in `examples` and parametrizes the `testcase` fixture for each discovered subdirectory.
    Ensures testcase scope is session.
    """
    return request.param


@pytest.fixture(scope="session")
def testcaserundir(rundir: Path, testgroup: str, testcase: str):
    """
    Generate the run directory at `{rundir}/{testgroup}/{testcase}` based on `examples/{testcase}`
    """
    # assert that `{rundir}` and `{rundir}/{testgroup}` exist
    if not rundir.exists():
        rundir.mkdir()
    if not (rundir / "data").exists():
        (rundir / "data").symlink_to(Path(__file__).parent / "data")
    if not (rundir / testgroup).exists():
        (rundir / testgroup).mkdir()
    # create the testcase directory
    sourcedir = Path(__file__).parent / "examples" / testcase
    targetdir = rundir / testgroup / testcase
    if targetdir.exists():
        shutil.rmtree(targetdir)
    shutil.copytree(sourcedir, targetdir, symlinks=True)
    # special treatment
    match testgroup:
        case "shortrun":
            helpers.adapt_parameter_file(
                sourcedir / "parameter.ini", targetdir / "parameter.ini", testlevel=2
            )
    return targetdir


# === TESTS === #


def test_examples(binpath, testgroup, testcaserundir, testcase, dryrun):
    """
    Test end2end GVEC runs with `{testgroup}/{testcase}/parameter.ini`

    Note: the `testgroup` fixture is contained in `testcaserundir`, but is given additionally to the test function
    for better readability and proper parameter ordering for the pytest nodeID.
    """
    args = [binpath, "parameter.ini"]
    # find restart statefile
    if "_restart" in testcase:
        base_name, suffix = re.match(
            r"(.*)(_restart.*)", testcase
        ).groups()  # remove sufix
        base_directory = testcaserundir / ".." / base_name
        states = [
            sd
            for sd in os.listdir(base_directory)
            if "State" in sd and sd.endswith(".dat")
        ]
        if not dryrun:
            args.append(base_directory / sorted(states)[-1])
        else:
            args.append(base_directory / "STATEFILE?")
    # run gvec
    with helpers.chdir(testcaserundir):
        if dryrun:
            with open("dryrun-examples", "w") as file:
                file.write(f"Dryrun would execute:\n{args}")
            return
        # run GVEC
        with open("stdout", "w") as stdout, open("stderr", "w") as stderr:
            subprocess.run(args, text=True, stdout=stdout, stderr=stderr)
        # check if GVEC was successful
        helpers.assert_empty_stderr()
        helpers.assert_stdout_finished()


@pytest.mark.regression
def test_regression(testgroup, testcase, rundir, refdir, dryrun):
    """
    Regression test of example GVEC runs and restarts.

    Compares all statefiles in `{rundir}/{testgroup}/{testcase}` with the corresponding statefiles in
    `{refdir}/{testgroup}/{testcase}` using `helpers.assert_equal_statefiles`.
    """
    directory = rundir / testgroup / testcase
    refdirectory = refdir / testgroup / testcase
    assert refdirectory.exists(), f"Reference testcase {refdirectory} does not exist"
    if dryrun:
        with open(directory / "dryrun-regression", "w") as file:
            file.write(
                f"Dryrun would compare directories:\n{directory}\nwith\n{refdirectory}\n"
            )
    # compare output to reference
    for filename in os.listdir(directory):
        if re.match(r"\w+State\w+\.dat", filename) and filename in os.listdir(
            refdirectory
        ):
            if not dryrun:
                helpers.assert_equal_statefiles(
                    directory / filename,
                    refdirectory / filename,
                )
            else:
                with open(directory / "dryrun-regression", "a") as file:
                    file.write(
                        f"Dryrun would compare files:\n{directory / filename}\nwith\n{refdirectory / filename}\n"
                    )
