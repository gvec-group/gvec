import sys
import os
import pytest
import logging
from pathlib import Path

# add helpers.py to the `pythonpath` to be importable by all tests
sys.path.append(str((Path(__file__).parent / "helpers")))


def pytest_addoption(parser):
    parser.addoption(
        "--builddir",
        type=Path,
        default=Path(__file__).parent / "../build",
        help="Path to builddir from git repo root folder",
    )
    parser.addoption(
        "--refdir",
        type=Path,
        default=None,
        help="Path to reference pytest root directory. Required for regression tests.",
    )


def pytest_configure(config):
    # register an additional marker
    config.addinivalue_line(
        "markers", "example: mark test as a example"
    )
    config.addinivalue_line(
        "markers", "shortrun: mark test as a shortrun"
    )
    config.addinivalue_line(
        "markers", "regression: mark test as a regression"
    )


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
    if request.config.getoption("--refdir") is None:
        pytest.exit("--refdir is required for regression tests")
    return Path(request.config.getoption("--refdir")).absolute()


@pytest.fixture(scope="session")
def rootdir(request) -> Path:
    """path to the (pytest root) directory"""
    return Path(request.config.rootdir).absolute()



@pytest.fixture
def logger(caplog):
    """get the pytest logger (with level set to DEBUG)"""
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger()
    yield logger
