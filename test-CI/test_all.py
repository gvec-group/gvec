import pytest
import os
import subprocess
import re
from pathlib import Path
import shutil

import helpers

try:
    import xarray as xr
    import gvec
except ImportError:
    pass


# === HELPER FUNCTIONS & FIXTURES === #


def discover_subdirs(path: Path | str = ".") -> list[str]:
    """
    Discover subdirectories in this module's directory for dynamic test detection.

    Ignores directories starting with `_` or `.`.

    Args:
        path (optional): Path to the directory to discover subdirectories in.\
            Defaults to the directory of this file.\
            Relative paths are relative to the directory of this file.

    Returns:
        list: A sorted list of subdirectory names.
    """
    # append relative path to the path of this file
    path = Path(path)
    if not path.is_absolute():
        path = Path(__file__).parent / path
    assert path.is_dir()
    # list subdirectories except for those starting with `_` or `.`
    return sorted(
        [
            sd
            for sd in os.listdir(path)
            if (path / sd).is_dir()
            and not sd.startswith("_")
            and not sd.startswith(".")
        ]
    )


@pytest.fixture(scope="session", params=discover_subdirs("examples"))
def testcase(request):
    """
    Fixture to parametrize the `testcase`

    Discoveres subdirectories in `examples` and parametrizes the `testcase` fixture for each discovered subdirectory.
    Ensures testcase scope is session.
    """
    return request.param


@pytest.fixture(scope="session", params=["example", "shortrun", "debugrun"])
def testgroup(request) -> str:
    """available test group names, will be automatically marked"""
    return request.param


# === TESTS === #


@pytest.mark.run_stage
def test_run(
    runargs_prefix,
    binpath,
    testgroup,
    testcaserundir,
    testcase,
    dryrun,
    annotations,
    artifact_pages_path,
):
    """
    Test end2end GVEC runs with `{testgroup}/{testcase}/parameter.ini`

    Note: the `testgroup` fixture is contained in `testcaserundir`, but is given additionally to the test function
    for better readability and proper parameter ordering for the pytest nodeID.
    """
    args = runargs_prefix + [binpath / "gvec", "parameter.ini"]
    # find restart statefile
    if "_restart" in testcase:
        base_name, suffix = re.match(
            r"(.*)(_restart.*)", testcase
        ).groups()  # remove sufix
        base_directory = testcaserundir / ".." / base_name
        if base_directory.exists():
            states = [
                sd
                for sd in os.listdir(base_directory)
                if "State" in sd and sd.endswith(".dat")
            ]
        if not dryrun:
            assert base_directory.exists(), "no base directory found for restart"
            assert (
                len(states) > 0
            ), f"no statefile for restart found in base directory ../{base_name}"
            laststatefile = sorted(states)[-1]
            (testcaserundir / laststatefile).symlink_to(
                os.path.relpath(base_directory / laststatefile, testcaserundir)
            )
            args.append(laststatefile)
        else:
            if base_directory.exists():
                if len(states) > 0:
                    laststatefile = sorted(states)[-1]
                else:
                    laststatefile = "STATEFILE???"
                args.append(laststatefile)
            else:
                args.append(base_directory / "???")
    # run gvec
    with helpers.chdir(testcaserundir):
        if dryrun:
            with open("dryrun-run.txt", "w") as file:
                file.write(f"DRYRUN: execute:\n {args} \n")
            return
        # run GVEC
        with open("stdout.txt", "w") as stdout:
            stdout.write(f"RUNNING: \n {args} \n")
        with open("stdout.txt", "a") as stdout, open("stderr.txt", "w") as stderr:
            subprocess.run(args, text=True, stdout=stdout, stderr=stderr)
        for filename in ["stdout", "stderr"]:
            if pages_rundir := os.environ.get("CASENAME"):
                pages_rundir = f"CIrun_{pages_rundir}"
            else:
                pages_rundir = "."
            annotations["gvec-output"].append(
                dict(
                    external_link=dict(
                        label=f"{testgroup}/{testcase}/{filename}",
                        url=f"{artifact_pages_path}/{pages_rundir}/{testgroup}/{testcase}/{filename}.txt",
                    )
                )
            )
        # check if GVEC was successful
        helpers.assert_empty_stderr(slurm="srun" in runargs_prefix)
        helpers.assert_stdout_finished(message="GVEC SUCESSFULLY FINISHED!")
        helpers.assert_stdout_OpenMP_MPI()


class BaseTestPost:
    """Base class for postprocessing/converter tests

    required fixtures/parameters:
    - name: identifier used in the directory path
    - args: list of arguments to pass to the executable
    """

    ciprefix = "CIpost"
    exec: str  # name of the executable

    @pytest.fixture(autouse=True)
    def testcasedir(
        self, postdir: Path, rundir: Path, testgroup: str, testcase: str, name: str
    ):
        """
        Generate the post directory at `{convdir}/{name}/{testgroup}/{testcase}` based on `{rundir}/{testgroup}/{testcase}`
        where name can include the identifier for different arguments
        """
        if not postdir.exists():
            postdir.mkdir()
        if not (postdir / name).exists():
            (postdir / name).mkdir()
        if not (postdir / name / "data").exists():
            (postdir / name / "data").symlink_to(Path(__file__).parent / "data")
        if not (postdir / name / testgroup).exists():
            (postdir / name / testgroup).mkdir()
        # create the testcase directory
        sourcedir = Path(__file__).parent / "examples" / testcase
        sourcerundir = rundir / testgroup / testcase
        targetdir = postdir / name / testgroup / testcase
        if targetdir.exists():
            shutil.rmtree(targetdir)
        # copy input files from examples/testcase
        shutil.copytree(sourcedir, targetdir, symlinks=True)
        # link to statefiles from run_stage
        states = sorted(sourcerundir.glob("*State*.dat"))
        for statefile in states:
            (targetdir / statefile.name).symlink_to(
                os.path.relpath(statefile, targetdir)
            )
        with helpers.chdir(targetdir):
            yield

    def test_post(
        self,
        testgroup,
        testcase,
        binpath,
        runargs_prefix,
        name,
        args,
        dryrun,
        annotations,
        artifact_pages_path,
    ):
        args = runargs_prefix + [binpath / self.exec] + args
        states = sorted(Path(".").glob("*State*.dat"))
        # dryrun
        if dryrun:
            with open("dryrun.txt", "w") as file:
                file.write(f"found statefiles {states}\n")
                file.write(f"DRYRUN: execute:\n {args} \n")
            return
        # insert the last statefile
        assert len(states) > 0, "no statefile for post/converter found"
        args[args.index("STATEFILE")] = states[-1]
        # run
        with open("stdout.txt", "w") as stdout:
            stdout.write(f"RUNNING: \n {args} \n")
        with (
            open("stdout.txt", "a") as stdout,
            open("stderr.txt", "w") as stderr,
        ):
            subprocess.run(args, text=True, stdout=stdout, stderr=stderr)
        # add link to artifact (CI)
        for filename in ["stdout", "stderr"]:
            if casename := os.environ.get("CASENAME"):
                pages_convdir = f"{self.ciprefix}_{casename}"
            else:
                pages_convdir = "."
            annotations["gvec-output"].append(
                dict(
                    external_link=dict(
                        label=f"{name}/{testgroup}/{testcase}/{filename}",
                        url=f"{artifact_pages_path}/{pages_convdir}/{name}/{testgroup}/{testcase}/{filename}.txt",
                    )
                )
            )


@pytest.mark.post_stage
class TestPost(BaseTestPost):
    exec = "gvec_post"

    @pytest.fixture(autouse=True)
    def adapt_parameters(self, util, testcasedir):
        # uncomment visualization flags
        util.adapt_parameter_file(
            "parameter.ini",
            "parameter.ini",
            visu1D="!0",
            visu2D="!0",
            visu3D="!0",
            SFLout="!-1",
        )

    @pytest.fixture
    def name(self):
        return "post"

    @pytest.fixture
    def args(self):
        return ["parameter.ini", "STATEFILE"]

    def test_post(
        self,
        testgroup,
        testcase,
        binpath,
        runargs_prefix,
        name,
        args,
        dryrun,
        annotations,
        artifact_pages_path,
    ):
        super().test_post(
            testgroup,
            testcase,
            binpath,
            runargs_prefix,
            name,
            args,
            dryrun,
            annotations,
            artifact_pages_path,
        )
        if dryrun:
            return

        helpers.assert_empty_stderr(slurm="srun" in runargs_prefix)
        helpers.assert_stdout_finished(message="GVEC POST FINISHED !")
        helpers.assert_stdout_OpenMP_MPI()


@pytest.mark.converter_stage
@pytest.mark.parametrize(
    "name, exec, msg, args",
    [
        (
            "to_gene",
            "test_gvec_to_gene",
            "TEST GVEC TO GENE",
            ["STATEFILE"],
        ),
        (
            "to_hopr",
            "test_gvec_to_hopr",
            "TEST GVEC TO HOPR",
            ["STATEFILE"],
        ),
    ]
    + [
        (
            f"to_castor3d-{i}",
            "convert_gvec_to_castor3d",
            "CONVERT GVEC TO CASTOR3D",
            args + ["STATEFILE", "gvec2castor3d_sfl.dat"],
        )
        for i, args in enumerate(
            [
                [
                    "--rpoints=7",
                    "--polpoints=12",
                    "--torpoints=8",
                    "--sflcoord=0",
                ],
                [
                    "--rpoints=8",
                    "--polpoints=11",
                    "--torpoints=9",
                    "--sflcoord=1",
                    "--factorsfl=2",
                ],
                [
                    "--rpoints=9",
                    "--polpoints=10",
                    "--torpoints=10",
                    "--sflcoord=2",
                    "--factorsfl=2",
                ],
                [
                    "--rpoints=6",
                    "--polpoints=10",
                    "--torpoints=10",
                    "--sflcoord=2",
                    "--factorsfl=2",
                    "--booz_relambda=0",
                ],
            ]
        )
    ]
    + [
        (
            "to_jorek",
            "convert_gvec_to_jorek",
            "CONVERT GVEC TO JOREK",
            [
                "--rpoints=8",
                "--npfactor=1",
                "--polpoints=12",
                "STATEFILE",
                "gvec2jorek_out.dat",
            ],
        )
    ],
    ids=[
        "to_gene",
        "to_hopr",
        "to_castor3d-0",
        "to_castor3d-1",
        "to_castor3d-2",
        "to_castor3d-3",
        "to_jorek",
    ],
)
class TestConverters(BaseTestPost):
    ciprefix = "CIconv"

    @pytest.fixture(autouse=True)
    def set_exec(self, exec):
        self.exec = exec

    def test_post(
        self,
        testgroup,
        testcase,
        binpath,
        runargs_prefix,
        name,
        args,
        msg,
        dryrun,
        annotations,
        artifact_pages_path,
    ):
        super().test_post(
            testgroup,
            testcase,
            binpath,
            runargs_prefix,
            name,
            args,
            dryrun,
            annotations,
            artifact_pages_path,
        )
        if dryrun:
            return

        helpers.assert_empty_stderr("stderr.txt", slurm="srun" in runargs_prefix)
        helpers.assert_stdout_finished("stdout.txt", message=f"{msg} FINISHED!")


@pytest.mark.pygvec
@pytest.mark.converter_stage
class TestToCAS3D(BaseTestPost):
    ciprefix = "CIconv"
    exec = "gvec_to_cas3d"

    @pytest.fixture
    def name(self):
        return "to_cas3d"

    @pytest.fixture
    def args(self):
        return [
            "--ns",
            "3",
            "--MN_out",
            "4",
            "4",
            "--stellsym",
            "--pointwise",
            "Booz-CAS3D.nc",
            "parameter.ini",
            "STATEFILE",
            "BoozFT-CAS3D.nc",
        ]

    def test_post(
        self,
        testgroup,
        testcase,
        binpath,
        runargs_prefix,
        name,
        args,
        dryrun,
        annotations,
        artifact_pages_path,
    ):
        super().test_post(
            testgroup,
            testcase,
            binpath,
            runargs_prefix,
            name,
            args,
            dryrun,
            annotations,
            artifact_pages_path,
        )
        if dryrun:
            return

        with open("stderr.txt") as file:
            lines = file.readlines()
            assert "done" in lines[-2]

        assert Path("BoozFT-CAS3D.nc").exists()
        boozft = xr.open_dataset("BoozFT-CAS3D.nc")
        for var in boozft.variables:
            assert "long_name" in boozft[var].attrs
            assert "symbol" in boozft[var].attrs

        assert Path("Booz-CAS3D.nc").exists()
        booz = xr.open_dataset("Booz-CAS3D.nc")
        for var in booz.variables:
            assert "long_name" in booz[var].attrs
            assert "symbol" in booz[var].attrs


@pytest.mark.regression_stage
def test_regression(
    testgroup,
    testcase,
    rundir,
    refdir,
    dryrun,
    logger,
    reg_rtol,
    reg_atol,
    extra_ignore_patterns,
):
    """
    Regression test of example GVEC runs and restarts.

    Compares all statefiles in `{rundir}/{testgroup}/{testcase}` with the corresponding statefiles in
    `{refdir}/{testgroup}/{testcase}` using `helpers.assert_equal_files`.
    """
    # the two directories to compare
    testcaserundir = rundir / testgroup / testcase
    testcaserefdir = refdir / testgroup / testcase
    # skip if any of the two directories do not exist
    if not testcaserundir.exists():
        logger.error("Testcase does not exist")
        pytest.fail("Testcase does not exist")
    if not testcaserefdir.exists():
        logger.error("Reference does not exist")
        pytest.fail("Reference does not exist")
    # compare the list of files in the two directories
    runfiles, reffiles = (
        set(
            [
                file
                for file in os.listdir(directory)
                if "dryrun" not in file and not file.endswith("~")
            ]
        )
        for directory in [testcaserundir, testcaserefdir]
    )
    # dry-run
    if dryrun:
        with open(testcaserundir / "dryrun-regression", "w") as file:
            file.write(
                f"DRYRUN: compare directories:\n{testcaserundir}\nwith\n{testcaserefdir}\n"
            )
            file.write(f"DRYRUN: compare files: {runfiles}\n")
        return
    # compare output to reference
    results = {}
    num_statefiles = 0
    num_diff_files = 0
    num_diff_lines = 0
    num_warnings = 0
    for filename in runfiles & reffiles:
        # skip symbolic links (mostly unreachable input)
        if (testcaserundir / filename).is_symlink():
            results[filename] = "symlink"
            continue
        # statefiles
        elif re.match(r".*_State_[\d_]*\.dat", filename):
            num = helpers.check_diff_files(
                testcaserundir / filename,
                testcaserefdir / filename,
                ignore_regexs=[r".*/.*"]
                + extra_ignore_patterns,  # ignore lines with a path
                atol=reg_atol,
                rtol=reg_rtol,
            )
            num_statefiles += 1
        elif re.match(r"log.*\.csv", filename):
            num = helpers.check_diff_files(
                testcaserundir / filename,
                testcaserefdir / filename,
                ignore_regexs=extra_ignore_patterns,
                ignore_columns=[1],
                atol=reg_atol,
                rtol=reg_rtol,
            )
        elif filename == "stdout.txt":
            num = helpers.check_diff_files(
                testcaserundir / filename,
                testcaserefdir / filename,
                ignore_regexs=[
                    r".* |- .*",
                    r".*GIT_.*",
                    r".*CMAKE.*",
                    r".*sec.*",
                    r".*date.*",
                    r".*PosixPath.*",
                    r"^[\s=]*$",
                    r"100%\| \.\.\. of",
                ]
                + extra_ignore_patterns,
                warn_regexs=[
                    "Number of OpenMP threads",
                    "Number of MPI tasks",
                    "GIT_",
                    "CMAKE",
                ],
                atol=reg_atol,
                rtol=reg_rtol,
            )
        else:
            results[filename] = "ignored"
            continue
        num_warnings += num[2]
        if num[0] > 0 or num[1] > 0:
            num_diff_files += 1
            num_diff_lines += num[0] + num[1]
            if num[0] and num[1]:
                results[filename] = "t&ndiff"
            elif num[0]:
                results[filename] = "txtdiff"
            else:
                results[filename] = "numdiff"
        else:
            results[filename] = "success"
    if num_statefiles == 0:
        logger.warning("Did not compare any statefiles!?")
        pytest.raised_warnings = True
    if num_warnings > 0:
        logger.warning(f"Found {num_warnings} warnings!")
        pytest.raised_warnings = True
    logger.info(f"{' SUMMARY ':=^80}")
    red, green, reset = "\x1b[31;20m", "\x1b[32;20m", "\x1b[0m"
    for filename, result in sorted(results.items(), key=lambda x: (x[1], x[0])):
        if result == "success":
            logger.info(f"... {green}{result}{reset} : {filename}")
        elif result in ["ignored", "symlink"]:
            logger.debug(f"... {result} : {filename}")
        else:
            logger.error(f"... {red}{result}{reset} : {filename}")
    for filename in runfiles - reffiles:
        logger.error(f"... {red}extra{reset}   : {filename}")
    for filename in reffiles - runfiles:
        logger.error(f"... {red}missing{reset} : {filename}")
    if num_diff_files > 0 or runfiles != reffiles:
        msg = (
            f"Found {num_diff_files} different files with {num_diff_lines} different lines, "
            f"{len(runfiles - reffiles)} additional files and {len(reffiles - runfiles)} missing files."
        )
        logger.info(f"{' SUMMARY ':=^80}")
        logger.error(msg)
        raise AssertionError(msg)
