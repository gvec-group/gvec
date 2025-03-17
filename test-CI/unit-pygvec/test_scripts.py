"""test the gvec.scripts subpackage

Only very simple tests for now.
"""

import os
from pathlib import Path
import shutil
import subprocess
import re

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
