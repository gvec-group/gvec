import pytest
import sys, os
import subprocess
import re

import helpers

# === AUTOMATIC / SETUP FIXTURES === #


@pytest.fixture(autouse=True, scope="module")
def chdir():
    """change working directory to this module's directory"""
    with helpers.chdir(os.path.dirname(__file__)):
        yield  # ToDo: does this work as expected?


# === HELPER FUNCTIONS === #


def discover_subdirs():
    """discover subdirectories with names `test_*` for dynamic test detection"""
    with helpers.chdir(os.path.dirname(__file__)):
        subdirs = os.listdir()
        return sorted(
            [sd for sd in subdirs if os.path.isdir(sd) and sd.startswith("test_")]
        )


# === TESTS === #


@pytest.mark.skip  # old / tryout
def test_try1():
    path = "test_shortrun/"
    name = "GITLAB_RUN_State_0000_00001000.dat"
    name2 = "GITLAB_RUN_State_0000_00001000.dat"
    helpers.assert_equal_statefiles(f"{path}/REF_{name}", f"{path}/{name2}")


#@pytest.mark.skip  # old / tryout
def test_run_parameter1():
    cmd = "../../../build/bin/gvec"
    cmdref = "../../../build/bin/gvec"
    with helpers.chdir("test_shortrun"):
        with helpers.chdir("ref"):
            pass
        subprocess.run([cmd, "parameter.ini"], capture_output=True)
        helpers.assert_equal_statefiles(
            "GITLAB_RUN_State_0000_00001000.dat",
            "REF_GITLAB_RUN_State_0000_00001000.dat",
        )


@pytest.mark.parametrize("directory", discover_subdirs())
def test_run_parameter(cmd, directory):
    """test full GVEC run with `parameter.ini`

    Automatically detects directories at the same level as this file,
    which start with `test_*` and tries to run GVEC on the `parameter.ini`
    within each directory. All (state) files starting with `REF*.dat` are
    compared with the respective GVEC output file.
    """
    with helpers.chdir(directory):
        with open("stdout", "w") as stdout, open("stderr", "w") as stderr:
            subprocess.run(
                [cmd, "parameter.ini"], text=True, stdout=stdout, stderr=stderr
            )
        # check if GVEC was successful
        helpers.assert_empty_stderr()
        helpers.assert_stdout_finished()
        # compare output to reference
        for path in os.listdir():
            if match := re.match(r"(REF_(\w+\.dat))", path):
                helpers.assert_equal_statefiles(
                    match.groups()[0],
                    match.groups()[1],
                )
