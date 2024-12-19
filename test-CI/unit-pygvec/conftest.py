import subprocess
import os

import pytest

import helpers

# === CONFIGURATION === #


def pytest_configure(config):
    """
    Add custom markers to pytest, as if using the `pytest.ini` file.

    Args:
        config (Config): The pytest configuration object.
    """
    # register additional markers
    for marker in [
        "unit: mark test as a unittest - a test that is isolated and tests a specific functionality",
        "pygvec: mark test that are related to the python bindings",
    ]:
        config.addinivalue_line("markers", marker)


def pytest_collection_modifyitems(session, config, items):
    """
    Modify collected tests

    Add global markers to each test.

    Args:
        items (List[Item]): The collected pytest testitem objects.
    """
    for item in items:
        if "unit-pygvec" in item.nodeid:
            item.add_marker(pytest.mark.unit)
            item.add_marker(pytest.mark.pygvec)
            if config.getoption("--dry-run"):
                item.add_marker(pytest.mark.skip("Dry-run: skipping unit tests"))


# === FIXTURES === #


@pytest.fixture(scope="session")
def testgroup():
    return "unit-pygvec"


@pytest.fixture(scope="session")
def testcase():
    return "frenet_axisNB_N2-12"


@pytest.fixture(scope="session")
def testcase_run(testgroup, testcaserundir, testcase, annotations, artifact_pages_path):
    """run the default testcase and store the output"""
    # assume pip-build: the gvec executable is in the PATH
    args = ["gvec", "parameter.ini"]
    # run gvec - adapted from test_all.test_run
    with helpers.chdir(testcaserundir):
        # skip if stdout & stderr are present and ok
        try:
            helpers.assert_empty_stderr()
            helpers.assert_stdout_finished(message="GVEC SUCESSFULLY FINISHED!")
        except (FileNotFoundError, AssertionError):
            pass
        else:
            return
        # run GVEC
        with open("stdout.txt", "w") as stdout:
            stdout.write(f"RUNNING: \n {args} \n")
        with open("stdout.txt", "a") as stdout, open("stderr.txt", "w") as stderr:
            subprocess.run(args, text=True, stdout=stdout, stderr=stderr)
        for filename in ["stdout", "stderr"]:
            if pages_rundir := os.environ.get("CASENAME"):
                pages_rundir = f"CIrun_{pages_rundir}"
            else:
                pages_rundir = "."
            annotations["gvec-output"].append(
                dict(
                    external_link=dict(
                        label=f"{testgroup}/{testcase}/{filename}",
                        url=f"{artifact_pages_path}/{pages_rundir}/{testgroup}/{testcase}/{filename}.txt",
                    )
                )
            )
        # check if GVEC was successful
        helpers.assert_empty_stderr()
        helpers.assert_stdout_finished(message="GVEC SUCESSFULLY FINISHED!")
    return


@pytest.fixture()
def testfiles(tmpdir, testcaserundir, testcase_run):
    """prepare the testcase parameters"""
    paramfile = "parameter.ini"
    statefile = "GVEC_axisNB_N2-12-hi_iota07_State_0000_00000001.dat"
    yield (testcaserundir / paramfile, testcaserundir / statefile)
