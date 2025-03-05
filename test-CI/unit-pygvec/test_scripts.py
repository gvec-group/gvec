"""test the gvec.scripts subpackage

Only very simple tests for now.
"""

import os
from pathlib import Path
import shutil

import pytest

try:
    import gvec
except ImportError:
    pass  # tests will be skipped via the `check_import` fixture


# === FIXTURES === #


@pytest.fixture(autouse=True)
def prepare_testcaserundir(tmp_path):
    """Prepare the test case run directory"""
    testcase = "w7x"
    shutil.copytree(
        Path(__file__).parent / "../examples/" / testcase, tmp_path, dirs_exist_ok=True
    )
    source = os.getcwd()
    os.chdir(tmp_path)
    yield
    os.chdir(source)


# === TESTS === #


def test_zero_current():
    """
    Test the zero_current script
    """
    from gvec.scripts.zero_current import run_zero_current

    run_zero_current(
        template="parameter.ini",
        max_iteration=2,
        gvec_miniter=5,
        gvec_maxiter=10,
        gvec_nelems_min=2,
        gvec_nelems_max=2,
        iota_poly_degree=5,
        reinit_LA=False,
        progressbar=False,
        outputfile="GVEC_zero_current.nc",
    )

    for i in range(2 + 1):
        assert Path(f"{i:02d}").exists()
        assert Path(f"{i:02d}/parameter.ini").exists()
        assert Path(f"{i:02d}/stdout.txt").exists()
        assert len(list(Path(f"{i:02d}").glob("*State*.dat"))) == 2

    assert list(Path("00").glob("*State*0000.dat"))
    assert list(Path("00").glob("*State*0005.dat"))
    assert list(Path("02").glob("*State*0000.dat"))
    assert list(Path("02").glob("*State*0010.dat"))
    assert Path("GVEC_zero_current.nc").exists()


def test_zero_current_cmd():
    """
    Test the zero_current script via command line
    """
    from subprocess import run

    run(
        [
            "gvec_zero_current",
            "parameter.ini",
            "--pi",
            "3",
            "--gi",
            "5",
            "10",
            "--ge",
            "2",
            "2",
            "-q",
        ]
    )

    for i in range(2 + 1):
        assert Path(f"{i:02d}").exists()
        assert Path(f"{i:02d}/parameter.ini").exists()
        assert Path(f"{i:02d}/stdout.txt").exists()
        assert len(list(Path(f"{i:02d}").glob("*State*.dat"))) == 2

    assert list(Path("00").glob("*State*0000.dat"))
    assert list(Path("00").glob("*State*0005.dat"))
    assert list(Path("02").glob("*State*0000.dat"))
    assert list(Path("02").glob("*State*0010.dat"))
    assert Path("GVEC_zero_current.nc").exists()
