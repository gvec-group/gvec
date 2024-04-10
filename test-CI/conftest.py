import sys
import os
import pytest
import logging
from pathlib import Path
import json
import re

# add helpers.py to the `pythonpath` to be importable by all tests
sys.path.append(str((Path(__file__).parent / "helpers")))


# === PYTEST CONFIGURATION === #


def pytest_addoption(parser):
    """
    Add custom command line options to pytest.

    Args:
        parser (ArgumentParser): The pytest argument parser.
    """

    group = parser.getgroup("custom_directories")
    group.addoption(
        "--builddir",
        type=Path,
        default=Path(__file__).parent / "../build",
        help="Path to build directory",
    )
    group.addoption(
        "--refdir",
        type=Path,
        default=None,
        help="Path to reference pytest root directory. Required for regression tests.",
    )
    group.addoption(
        "--rundir",
        type=Path,
        default=Path(__file__).parent / "run",
        help="Path to run directory",
    )
    group.addoption(
        "--postdir",
        type=Path,
        default=Path(__file__).parent / "post",
        help="Path to post directory",
    )
    group.addoption(
        "--annotations",
        type=Path,
        default=None,
        help="Path to the (output) annotations report file",
    )

    group = parser.getgroup("custom_run_options")
    group.addoption(
        "--dry-run",
        action="store_true",
        help="Do not execute commands, but generate all runs.",
    )
    group.addoption(
        "--reg-rtol",
        type=float,
        default=1.0e-7,
        help="relative tolerance for regression stage",
    )
    group.addoption(
        "--reg-atol",
        type=float,
        default=1.0e-10,
        help="absolute tolerance for regression stage",
    )

    group = parser.getgroup("custom_regression_options")
    group.addoption(
        "--add-ignore-pattern",
        action="append",
        default=[],
        help="additional ignore patterns for the `test_regression` line comparison",
    )
    group = parser.getgroup("custom_exit_code")
    group.addoption(
        "--suppress-no-test-exit-code",
        action="store_true",
        default=True,
        help='Suppress the "no tests collected" exit code.',
    )
    group.addoption(
        "--suppress-tests-failed-exit-code",
        action="store_true",
        default=False,
        help='Suppress the "some tests failed" exit code.',
    )


def pytest_configure(config):
    """
    Add custom markers to pytest, as if using the `pytest.ini` file.

    Args:
        config (Config): The pytest configuration object.
    """
    # register additional markers
    for marker in [
        "example: mark test as a example, which are all tests specified in `test-CI/examples`, executed in folder `rundir/example`",
        "restart: mark test as a restart (deduced from example folder name). Needs the example to be run first!",
        "shortrun: mark test as a shortrun  (executed in folder `rundir/shortrun`, overwrites parameters: `testlevel=-1` and `MaxIter=1`)",
        "debugrun: mark test as a debugrun (executed in folder `rundir/debugrun`, overwrites parameters: `testlevel=2` and `MaxIter=1`)",
        "run_stage: mark test belonging to the run stage (executed for all testgroups into a `rundir`)",
        "post_stage: mark test belonging to the post-processing stage (executed for all testgroups into a `postdir`, activates visualization parameters). Needs run_stage to be executed before in a given `rundir` directory.",
        "regression_stage: mark test belonging to the regression stage (compares files from `rundir` and  `refdir`. The --refdir argument is mandatory!",
    ]:
        config.addinivalue_line("markers", marker)

    # custom global variables
    pytest.raised_warnings = False


def pytest_collection_modifyitems(items):
    """
    Modify collected tests

    Add markers to each test based on the value of the `testgroup` (`example` or `shortrun`  or `debugrun`).

    Args:
        items (List[Item]): The collected pytest testitem objects.
    """
    for item in items:
        if "testgroup" in getattr(item, "fixturenames", ()):
            item.add_marker(getattr(pytest.mark, item.callspec.getparam("testgroup")))
        if ("testcase" in getattr(item, "fixturenames", ())) and (
            "_restart" in item.callspec.getparam("testcase")
        ):
            item.add_marker(getattr(pytest.mark, "restart"))
    # sort tests by testgroup and testcase
    stages = ["test_run", "test_regression", "test_post"]
    def sort_key(item):
        if not hasattr(item, "callspec") or "testgroup" not in item.callspec.params:
            return -1, item.name
        return stages.index(item.name.split("[")[0]), item.callspec.getparam("testgroup"), item.callspec.getparam("testcase")
    items.sort(key=sort_key)


def pytest_runtest_setup(item):
    """
    Runs before each test is executed.

    Skip regression tests if `--refdir` is not set

    Args:
        item (Item): The pytest testitem object.
    """
    if len(list(item.iter_markers(name="regression"))) and not item.config.getoption(
        "--refdir"
    ):
        pytest.skip("regression tests require `--refdir`")


def pytest_sessionfinish(session, exitstatus):
    # Original code taken from "pytest-custom_exit_code" plugin
    # From pytest version >=5, the values are inside an enum
    from pytest import ExitCode

    # === Custom exit codes === #
    no_tests_collected = ExitCode.NO_TESTS_COLLECTED  # 5
    tests_failed = ExitCode.TESTS_FAILED  # 1
    ok = ExitCode.OK  # 0
    if session.config.getoption("--suppress-no-test-exit-code"):
        if exitstatus == no_tests_collected:
            print(
                f"EXITING OK INSTEAD OF 'no tests collected' (exit status=0 instead of {exitstatus})"
            )
            session.exitstatus = ok

    if session.config.getoption("--suppress-tests-failed-exit-code"):
        if exitstatus == tests_failed:
            print(
                f"EXITING OK INSTEAD OF 'tests failed' (exit status=0 instead of {exitstatus})"
            )
            session.exitstatus = ok

    if pytest.raised_warnings:
        session.exitstatus = no_tests_collected
    


# === FIXTURES === #


@pytest.fixture(scope="session")
def builddir(request) -> Path:
    """path to the build directory"""
    return Path(request.config.getoption("--builddir")).absolute()


@pytest.fixture(scope="session")
def binpath(builddir) -> Path:
    """path to the binary folder"""
    return builddir / "bin"


@pytest.fixture(scope="session")
def dryrun(request) -> bool:
    """flag for dry-run mode"""
    return request.config.getoption("--dry-run")


@pytest.fixture(scope="session")
def reg_rtol(request) -> float:
    """relative tolerance for regression"""
    return request.config.getoption("--reg-rtol")


@pytest.fixture(scope="session")
def reg_atol(request) -> float:
    """absolute tolerance for regression"""
    return request.config.getoption("--reg-atol")


@pytest.fixture(scope="session")
def refdir(request) -> Path:
    """path to the reference (test-CI) directory"""
    if request.config.getoption("--refdir") is None:
        pytest.skip("--refdir is required for regression tests")
    return Path(request.config.getoption("--refdir")).absolute()


@pytest.fixture(scope="session")
def postdir(request) -> Path:
    """path to the post directory, default is test-CI/post"""
    return Path(request.config.getoption("--postdir")).absolute()


@pytest.fixture(scope="session")
def rootdir(request) -> Path:
    """path to the root (test-CI) directory"""
    return Path(request.config.rootdir).absolute()


@pytest.fixture(scope="session")
def rundir(request) -> Path:
    """path to the run directory for the test results"""
    return Path(request.config.getoption("--rundir")).absolute()


@pytest.fixture(scope="session")
def extra_ignore_patterns(request) -> list:
    """additional ignore patterns for the `test_regression` line comparison"""
    return request.config.getoption("--add-ignore-pattern")


@pytest.fixture(scope="session", params=["example", "shortrun", "debugrun"])
def testgroup(request) -> str:
    """available test group names, will be automatically marked"""
    return request.param


@pytest.fixture(scope="session")
def artifact_pages_path(rundir) -> str:
    """path to the CI artifact pages"""
    env = os.environ
    try:
        project_path = "/".join(env["CI_PROJECT_PATH"].split("/")[1:])  # remove the root namespace of the project path
        return f"{env['CI_SERVER_PROTOCOL']}://{env['CI_PROJECT_ROOT_NAMESPACE']}.{env['CI_PAGES_DOMAIN']}/-/{project_path}/-/jobs/{env['CI_JOB_ID']}/artifacts"
    except KeyError:
        return rundir


@pytest.fixture(scope="session", autouse=True)
def annotations(request, artifact_pages_path):
    """annotations for the CI report"""
    mapping = {}
    # load the annotations from the file if it exists
    if (path := request.config.getoption("--annotations")) and os.path.exists(path):
        with open(path, "r") as file:
            mapping = json.load(file)
    # required collections
    for key in ["pytest-log", "gvec-output"]:
        if key not in mapping:
            mapping[key] = []
    # add the pytest log file
    if path := request.config.getoption("--log-file"):
        mapping["pytest-log"].append(dict(external_link=dict(
            label = path,
            url = f"{artifact_pages_path}/{path}",
        )))
    # hand over mapping
    yield mapping
    # save the annotations to the file
    if path := request.config.getoption("--annotations"):
        with open(path, "w") as file:
            json.dump(mapping, file)


@pytest.fixture(autouse=True)
def logger(caplog):
    """get the pytest logger (with level set to DEBUG)"""
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger()
    return logger
