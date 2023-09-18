# THIS FILE IS USED TO CONFIGURE PYTEST WITHIN THIS REPOSITORY

import pytest

def pytest_addoption(parser):
    parser.addoption('--builddir',type=str, default='build',
                                 help='Path to builddir from git repo root folder, mandatory last argument')
    parser.addoption('--mpiruncmd',type=str, default='mpirun -n',
                                 help='start command when running in MPI, followed by nMPItasks.')
    parser.addoption('--nmpitasks',type=int, default=0,
                                 help='0: (default) do not run with MPI. >0: number of mpi tasks.')


@pytest.fixture(scope="session")
def builddir_option(request):
    ''' returns the value of option "--builddir" '''
    return request.config.getoption("--builddir")


def pytest_configure(config):
    print("config.option: ", config.option)


def pytest_sessionstart(session):
    import subprocess
    # check if numdiff exists!
    try:
        sbp=subprocess.run(['numdiff','--version'],capture_output=True,text=True)
        if(sbp.returncode ==0):
            print('\n check: numdiff found \n')
        else:
            pytest.exit('\n'+sbp.stdout+sbp.stderr+'\n ==> Error in pytest sessionstart: "numdiff" command not found!! \n')
    except:
        pytest.exit('\n ==> Error in pytest sessionstart: "numdiff" command not found!! \n')
