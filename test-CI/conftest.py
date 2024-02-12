import sys
import pytest
import logging
from pathlib import Path




# add helpers.py to the `pythonpath` to be importable by all tests
sys.path.append(str((Path(__file__).parent / "helpers")))


def pytest_addoption(parser):
    """
    Add custom command line options to pytest.

    Args:
        parser (ArgumentParser): The pytest argument parser.
    """
    parser.addoption(
        "--builddir",
        type=Path,
        default=Path(__file__).parent / "../build",
        help="Path to build directory",
    )
    parser.addoption(
        "--refdir",
        type=Path,
        default=None,
        help="Path to reference pytest root directory. Required for regression tests.",
    )
    parser.addoption(
        "--rundir",
        type=Path,
        default=Path(__file__).parent / "run",
        help="Path to run directory",
    )
    parser.addoption(
        "--dry-run",
        action="store_true",
        help="Do not execute commands, but generate all runs.",
    )


def pytest_configure(config):
    """
    Add custom markers to pytest, as if using the `pytest.ini` file.

    Args:
        config (Config): The pytest configuration object.
    """
    # register an additional marker
    config.addinivalue_line(
        "markers", "example: mark test as a example"
    )
    config.addinivalue_line(
        "markers", "shortrun: mark test as a shortrun  (overwrites parameters: `testlevel=-1` and `MaxIter=1`)"
    )
    config.addinivalue_line(
        "markers", "debugrun: mark test as a debugrun (overwrites parameters: `testlevel=2` and `MaxIter=1`) "
    )
    config.addinivalue_line(
        "markers", "regression: mark test as a regression"
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
    # sort tests by testgroup and testcase
    items.sort(key=lambda item: (item.callspec.getparam("testgroup"), item.callspec.getparam("testcase")))


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
    """path to the gvec binary"""
    return builddir / "bin" / "gvec"


@pytest.fixture(scope="session")
def dryrun(request) -> bool:
    """flag for dry-run mode"""
    return request.config.getoption("--dry-run")


@pytest.fixture(scope="session")
def refdir(request) -> Path:
    """path to the reference (test-CI) directory"""
    if request.config.getoption("--refdir") is None:
        pytest.exit("--refdir is required for regression tests")
    return Path(request.config.getoption("--refdir")).absolute()


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


@pytest.fixture
def logger(caplog):
    """get the pytest logger (with level set to DEBUG)"""
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger()
    yield logger
