"""test the gvec.scripts subpackage

Only very simple tests for now.
"""

import os
from pathlib import Path
import shutil
import subprocess
import re

import pytest
import numpy as np

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


def test_version():
    """
    Test the version of the main script
    """
    proc = subprocess.run(["pygvec", "-V"], capture_output=True)
    m = re.match(r"pyGVEC v(\S+)", proc.stdout.decode())
    assert m is not None
    assert m.group(1) == gvec.__version__
    assert proc.returncode == 0


@pytest.mark.parametrize("mode", ["", "run", "to-cas3d", "convert-params"])
def test_help(mode):
    """
    Test the help message of the main script
    """
    if mode == "":
        args = ["pygvec", "-h"]
    else:
        args = ["pygvec", mode, "-h"]

    proc = subprocess.run(args, capture_output=True)
    assert proc.returncode == 0


@pytest.mark.parametrize("suffix", ["ini", "yaml", "toml"])
def test_run_stages(suffix):
    """
    Test the run_stages function
    """
    args = [f"parameter.{suffix}"]
    gvec.scripts.run.main(args)

    if suffix == "ini":
        assert Path("W7X_State_0000_00000100.dat").exists()
    else:
        assert Path("0-00/W7X_State_0000_00000010.dat").exists()
        assert Path("0-01/W7X_State_0001_00000010.dat").exists()
        assert Path("1-00/W7X_State_0002_00000005.dat").exists()


def test_quasr_real_dft():
    def exfunc(x):
        return x * 0 + 3 + 1.4 * np.sin(2 * x + 0.4) + 0.3 * np.cos(4 * x - 0.3)

    def exfuncd(x):
        return x * 0 + 2 * 1.4 * np.cos(2 * x + 0.4) - 4 * 0.3 * np.sin(4 * x - 0.3)

    def exfuncdd(x):
        return x * 0 - 4 * 1.4 * np.sin(2 * x + 0.4) - 16 * 0.3 * np.cos(4 * x - 0.3)

    nzeta_test = 9
    nzeta_up = 14

    zeta_test = np.linspace(
        0, np.pi, nzeta_test, endpoint=False
    )  # data on one field period nfp=2 -> modes must be multiples of 2...
    zeta_up = np.linspace(0, 2 * np.pi, nzeta_up, endpoint=False)

    f1 = exfunc(zeta_test)

    rdft = gvec.scripts.quasr.real_dft_mat(zeta_test, zeta_up, nfp=2)
    f3 = rdft["BF"] @ f1

    d_rdft = gvec.scripts.quasr.real_dft_mat(zeta_test, zeta_up, deriv=1, nfp=2)

    df3 = (d_rdft["B"] @ (d_rdft["F"] @ f1)).real

    dd_rdft = gvec.scripts.quasr.real_dft_mat(zeta_test, zeta_up, deriv=2, nfp=2)

    ddf3 = dd_rdft["BF"] @ f1

    assert np.allclose(f3, exfunc(zeta_up))
    assert np.allclose(df3, exfuncd(zeta_up))
    assert np.allclose(ddf3, exfuncdd(zeta_up))
