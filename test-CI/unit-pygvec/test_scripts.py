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
    import numpy as np
    import xarray as xr

    import gvec
except ImportError:
    pass  # tests will be skipped via the `check_import` fixture


DATA = Path(__file__).parent / "../data"
QUASR_ID = 112714

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

    match suffix:
        case "ini":
            assert Path("W7X_State_0000_00000100.dat").exists()
        case "yaml":
            assert Path("W7X_yaml_gvec_stages").exists()
            assert Path("W7X_yaml_State_final.dat").exists()
            assert Path("parameter_W7X_yaml_final.ini").exists()
        case "toml":
            assert Path("W7X_toml_gvec_stages").exists()
            assert Path("W7X_toml_State_final.dat").exists()
            assert Path("parameter_W7X_toml_final.ini").exists()


def test_picard_auto():
    """
    Test the if the picard_current auto mode achieves "zero" current
    """
    parameters = gvec.util.read_parameters("parameter.toml")
    parameters["picard_current"] = "auto"
    ProjectName = parameters["ProjectName"]
    if "stages" in parameters:
        with pytest.raises(ValueError):
            run_with_stages = gvec.scripts.run.RunWithStages(parameters)
        del parameters["stages"]
    run_with_stages = gvec.scripts.run.RunWithStages(parameters)
    rundir, final_state, diagnostics = run_with_stages.run_stages(parameters)
    assert diagnostics.force_X1[-1].data <= 1e-4
    assert diagnostics.force_X2[-1].data <= 1e-4
    assert diagnostics.force_LA[-1].data <= 1e-4

    assert Path(f"{ProjectName}_State_final.dat").exists()
    assert Path(f"parameter_{ProjectName}_final.ini").exists()
    rho = np.sqrt(np.linspace(0, 1, 101))
    rho[0] = 1e-4
    with gvec.State(f"parameter_{ProjectName}_final.ini", final_state) as state:
        ev = gvec.Evaluations(rho=rho, theta="int", zeta="int", state=state)
        state.compute(ev, "I_tor")
    assert ev.I_tor.max().data < 1e-6


def test_stages_without_current():
    """Test if stages run without current constraint"""
    parameters = gvec.util.read_parameters("parameter.toml")
    parameters["picard_current"] = "off"
    ProjectName = parameters["ProjectName"]
    del parameters["I_tor"]
    parameters["stages"] = [
        {"minimize_tol": 1e-2, "sgrid": {"nelems": 2}},
        {"minimize_tol": 1e-3, "sgrid": {"nelems": 3}},
    ]
    run_with_stages = gvec.scripts.run.RunWithStages(parameters)
    rundir, final_state, diagnostics = run_with_stages.run_stages(parameters)
    assert Path(f"{ProjectName}_State_final.dat").exists()
    assert Path(f"parameter_{ProjectName}_final.ini").exists()


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


def test_quasr_download(tmp_path):
    json = gvec.scripts.quasr.get_json_from_quasr(
        QUASR_ID, tmp_path / "quasr-{QUASR_ID:07d}.json"
    )
    json_ref = DATA / f"quasr-{QUASR_ID:07d}.json"
    with (
        open(json, "r") as file,
        open(json_ref, "r") as reference,
    ):
        assert file.read().strip() == reference.read().strip()


def test_quasr_noerror(tmp_path):
    """
    Test the load-quasr script
    """
    hmap = Path(f"quasr-{QUASR_ID:07d}-Gframe.nc")
    with gvec.util.chdir(tmp_path):
        args = ["-f", str(DATA / f"quasr-{QUASR_ID:07d}-boundary.nc")]
        gvec.scripts.quasr.main(args)
        assert hmap.exists()
