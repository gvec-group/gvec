import sys
import os
import pytest
import logging
from pathlib import Path

# add helpers.py to the `pythonpath` to be importable by all tests
sys.path.append(os.path.join(os.path.dirname(__file__), "helpers"))


def pytest_addoption(parser):
    parser.addoption(
        "--builddir",
        type=str,
        default="build",
        help="Path to builddir from git repo root folder",
    )
    parser.addoption(
        "--refdir",
        type=str,
        default=None,
        help="Path to reference pytest root directory. Required for regression tests.",
    )
    parser.addoption(
        "-S", 
        "--stage",
        action="store",
        choices=["examples", "regression"],
        default=[],
        help="Test stage to run",
    )


def pytest_configure(config):
    # register an additional marker
    config.addinivalue_line(
        "markers", "example: mark test as a testable example"
    )


def pytest_runtest_setup(item):  # ToDo document
    if len(list(item.iter_markers("example"))) and "examples" not in item.config.getoption("-S"):
        pytest.skip("skipping examples test")
    elif "regression" in item.name and "regression" not in item.config.getoption("-S"):
        pytest.skip("skipping regression test")


@pytest.fixture(scope="session")
def builddir(request) -> Path:
    """path to the build directory"""
    return Path(request.config.getoption("--builddir")).absolute()


@pytest.fixture(scope="session")
def binpath(builddir) -> Path:
    """path to the gvec binary"""
    return builddir / "bin" / "gvec"


@pytest.fixture(scope="session")
def refdir(request) -> Path:
    """path to the reference (pytest root) directory"""
    return Path(request.config.getoption("--refdir")).absolute()


@pytest.fixture(scope="session")
def rootdir(config) -> Path:
    """path to the (pytest root) directory"""
    return Path(config.rootdir).absolute()



@pytest.fixture
def logger(caplog):
    """get the pytest logger (with level set to DEBUG)"""
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger()
    yield logger
