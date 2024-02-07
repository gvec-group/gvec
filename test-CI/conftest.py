import sys
import os
import pytest
import logging

# add helpers.py to the `pythonpath` to be importable by all tests
sys.path.append(os.path.join(os.path.dirname(__file__), "helpers"))


def pytest_addoption(parser):
    parser.addoption(
        "--cmd",
        type=os.path.abspath,
        default="build/bin/gvec",
        help="path to the GVEC executable",
    )

@pytest.fixture(scope="session")
def cmd(request):
    """path to the GVEC executable"""
    return request.config.getoption("--cmd")


@pytest.fixture
def logger(caplog):
    """get the pytest logger (with level set to DEBUG)"""
    caplog.set_level(logging.DEBUG)
    logger = logging.getLogger()
    yield logger
