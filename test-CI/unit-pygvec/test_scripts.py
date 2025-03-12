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


@pytest.mark.parametrize("suffix", ["ini", "yaml", "toml"])
def test_zero_current_cmd(suffix):
    """
    Test the zero_current script via command line
    """
    from subprocess import run

    proc = run(
        [
            "pygvec",
            "run",
            f"parameter.{suffix}",
        ],
    )
    assert proc.returncode == 0

    if suffix == "ini":
        assert Path("W7X_State_0000_00000100.dat").exists()
    else:
        assert Path("0-00/W7X_State_0000_00000010.dat").exists()
        assert Path("0-01/W7X_State_0001_00000010.dat").exists()
        assert Path("1-00/W7X_State_0002_00000005.dat").exists()
