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
def prepare_testcaserundir(tmp_path, util):
    """Prepare the test case run directory"""
    testcase = "w7x"
    shutil.copytree(
        Path(__file__).parent / "../examples/" / testcase, tmp_path, dirs_exist_ok=True
    )
    source = Path.cwd()
    with util.chdir(tmp_path):
        yield


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

    assert Path(f"W7X_{suffix}_State_final.dat").exists()
    assert Path(f"parameter_W7X_{suffix}_final.ini").exists()


def test_picard_auto():
    """
    Test the if the picard_current auto mode achieves "zero" current
    """
    parameters = gvec.util.read_parameters("parameter.toml")
    parameters["picard_current"] = "auto"
    ProjectName = "Test_auto"
    parameters["projectname"] = ProjectName
    if "stages" in parameters:
        with pytest.raises(ValueError):
            run_with_stages = gvec.run(parameters)
        del parameters["stages"]
    run_with_stages = gvec.run(parameters)
    diagnostics = run_with_stages.diagnostics_minimizer
    assert diagnostics.force_X1[-1].data <= 1e-4
    assert diagnostics.force_X2[-1].data <= 1e-4
    assert diagnostics.force_LA[-1].data <= 1e-4

    assert Path(f"{ProjectName}_State_final.dat").exists()
    assert Path(f"parameter_{ProjectName}_final.ini").exists()
    rho = np.sqrt(np.linspace(0, 1, 101))
    rho[0] = 1e-4
    ev = run_with_stages.state.evaluate("I_tor", "iota_curr", rho=rho, theta="int", zeta="int")
    assert ev.iota_curr.max().data < 1e-6


@pytest.mark.parametrize("ptype", ["interpolation", "polynomial", "bspline"])
def test_I_tor_types(ptype):
    """Test if all types of I_tor profiles result in valid current constraints."""
    parameters = gvec.util.read_parameters("parameter.toml")
    parameters["picard_current"] = "auto"
    ProjectName = f"Test_Itor_type_{ptype}"
    parameters["projectname"] = ProjectName

    # integrated two power profile (1-x²)
    coefs = 6000 * np.array([0, 1, 0, 1 / 3])

    if "stages" in parameters:
        del parameters["stages"]

    match ptype:
        case "interpolation":
            rho2 = np.linspace(1e-4, 1, 50)
            # integrated two power profile (1-x²)
            I_tor = 6000 * (rho2 - rho2**3 / 3)
            parameters["I_tor"] = dict(type=ptype, rho2=rho2, vals=I_tor)
        case "polynomial":
            parameters["I_tor"] = dict(type=ptype, coefs=coefs)
        case "bspline":
            from .test_profiles import poly2bspl_coeff

            c_bspl = np.zeros(len(coefs))
            knots = np.concatenate([np.zeros(len(coefs)), np.ones(len(coefs))])
            for j in range(len(coefs)):
                c_bspl[j] = poly2bspl_coeff(coefs, j, knots)
            parameters["I_tor"] = dict(type=ptype, coefs=c_bspl, knots=knots)

    run_with_stages = gvec.run(parameters)
    diagnostics = run_with_stages.diagnostics_minimizer
    assert diagnostics.force_X1[-1].data <= 1e-4
    assert diagnostics.force_X2[-1].data <= 1e-4
    assert diagnostics.force_LA[-1].data <= 1e-4

    assert Path(f"{ProjectName}_State_final.dat").exists()
    assert Path(f"parameter_{ProjectName}_final.ini").exists()

    I_tor_rms = np.sqrt(
        (run_with_stages.diagnostics_run.I_tor_delta.isel(run=-1) ** 2).mean(dim="rad")
    )
    assert I_tor_rms.data < 1e-6, f"Expected ΔI_tor < 1e-6, got ΔI_tor:{I_tor_rms.data}"


def test_maxRestarts():
    """Test if maxRestarts aborts correctly"""
    parameters = gvec.util.read_parameters("parameter.toml")
    parameters["picard_current"] = gvec.util.CaseInsensitiveDict(
        dict(maxRestarts=10, iota_tol=1e-6, target="iota")
    )
    ProjectName = "Test_maxRestarts"
    parameters["projectname"] = ProjectName

    parameters["stages"] = [
        {"picard_current": {"maxRestarts": 0}},
        {
            "picard_current": {"maxRestarts": 1, "target": "iota_and_force"},
            "minimize_tol": 1e-4,
        },
        {"picard_current": {"maxRestarts": 2}},
    ]
    run_with_stages = gvec.run(parameters)

    for n, runs in enumerate(run_with_stages.n_runs_in_stage):
        assert n >= runs - 1, (
            f"In stage {n} maxRestarts was violated! allowed restarts: {n + 1}, performed restarts: {runs}"
        )


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
    run_with_stages = gvec.run(parameters)
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


def test_quasr_noerror(tmp_path, util):
    """
    Test the load-quasr script
    """
    hmap = Path(f"quasr-{QUASR_ID:07d}-Gframe.nc")
    with util.chdir(tmp_path):
        args = ["-f", str(DATA / f"quasr-{QUASR_ID:07d}-boundary.nc")]
        gvec.scripts.quasr.main(args)
        assert hmap.exists()
