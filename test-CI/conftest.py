import pytest

def pytest_addoption(parser):
    parser.addoption('--builddir',type=str, default='build',
                                 help='Path to builddir from git repo root folder, mandatory last argument')


@pytest.fixture(scope="session")
def builddir_option(request):
    ''' returns the value of option "--builddir" '''
    return request.config.getoption("--builddir")

def pytest_configure(config):
    print("config.option: ", config.option)