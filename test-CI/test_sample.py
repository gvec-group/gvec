import pytest
import helpers

def func(x):
    return x + 1


def xxtest_wrong_answer():
    assert func(3) == 5


@pytest.mark.parametrize('N', [100,4,5,7])
def test_realanswer(N):
    assert func(N) == N+1


def test_print_b(builddir_option):
    print ("Displaying builddir: %s" % builddir_option)
