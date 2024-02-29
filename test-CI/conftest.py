import sys
import pytest
import logging
from pathlib import Path

#### for VScode debugging
#import debugpy
#debugpy.listen(5678)
#debugpy.wait_for_client()
####


# add helpers.py to the `pythonpath` to be importable by all tests
sys.path.append(str((Path(__file__).parent / "helpers")))


def pytest_addoption(parser):
    """
    Add custom command line options to pytest.

    Args:
        parser (ArgumentParser): The pytest argument parser.
    """
    
    group = parser.getgroup('custom_directories')
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
    group = parser.getgroup('custom_run_options')
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
    group = parser.getgroup('custom_exit_code')
    group.addoption(
        '--suppress-no-test-exit-code',
        action='store_true',
        default=True,
        help='Suppress the "no tests collected" exit code.'
    )
    group.addoption(
        '--suppress-tests-failed-exit-code',
        action='store_true',
        default=False,
        help='Suppress the "some tests failed" exit code.'
    )



def pytest_configure(config):
    """
    Add custom markers to pytest, as if using the `pytest.ini` file.

    Args:
        config (Config): The pytest configuration object.
    """
    # register an additional marker
    config.addinivalue_line(
        "markers", "example: mark test as a example, which are all tests specified in `test-CI/examples`, executed in folder `rundir/example`"
    )
    config.addinivalue_line(
        "markers", "restart: mark test as a restart (deduced from example folder name). Needs the example to be run first!"
    )
    config.addinivalue_line(
        "markers", "shortrun: mark test as a shortrun  (executed in folder `rundir/shortrun`, overwrites parameters: `testlevel=-1` and `MaxIter=1`)"
    )
    config.addinivalue_line(
        "markers", "debugrun: mark test as a debugrun (executed in folder `rundir/debugrun`, overwrites parameters: `testlevel=2` and `MaxIter=1`) "
    )
    config.addinivalue_line(
        "markers", "run_stage: mark test belonging to the run stage (executed for all testgroups into a `rundir`)"    
    )

    config.addinivalue_line(
        "markers", "post_stage: mark test belonging to the post-processing stage (executed for all testgroups into a `postdir`, activates visualization parameters). Needs run_stage to be executed before in a given `rundir` directory. "    
    )
    config.addinivalue_line(
        "markers", "regression_stage: mark test belonging to the regression stage (compares files from `rundir` and  `refdir`. The --refdir argument is mandatory!"    
    )


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
        if ("testcase" in getattr(item, "fixturenames", ())) and ("_restart" in item.callspec.getparam("testcase")):
            item.add_marker(getattr(pytest.mark, "restart"))
    # sort tests by testgroup and testcase
    stages = ["test_run", "test_regression", "test_post"]
    items.sort(key=lambda item: (stages.index(item.name.split('[')[0]), item.callspec.getparam("testgroup"), item.callspec.getparam("testcase")))


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


@pytest.fixture(scope="session")
def builddir(request) -> Path:
    """path to the build directory"""
    return Path(request.config.getoption("--builddir")).absolute()


@pytest.fixture(scope="session")
def binpath(request, builddir) -> Path:
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
def refdir(request,dryrun) -> Path:
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


@pytest.fixture(scope="session", params=["example", "shortrun","debugrun"])
def testgroup(request) -> str:
    """available test group names, will be automatically marked"""
    return request.param


@pytest.fixture(autouse=True)
def logger(caplog):
    """get the pytest logger (with level set to DEBUG)"""
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger()
    return logger

@pytest.hookimpl()
def pytest_sessionfinish(session, exitstatus):
    # Original code taken from "pytest-custom_exit_code" plugin
    # From pytest version >=5, the values are inside an enum
    from pytest import ExitCode
    no_tests_collected = ExitCode.NO_TESTS_COLLECTED
    tests_failed = ExitCode.TESTS_FAILED
    ok = ExitCode.OK
    if session.config.getoption('--suppress-no-test-exit-code'):
        if exitstatus == no_tests_collected:
            print(f"EXITING OK INSTEAD OF 'no tests collected' (exit status=0 instead of {exitstatus})")
            session.exitstatus = ok

    if session.config.getoption('--suppress-tests-failed-exit-code'):
        if exitstatus == tests_failed:
            print(f"EXITING OK INSTEAD OF 'tests failed' (exit status=0 instead of {exitstatus})")
            session.exitstatus = ok
