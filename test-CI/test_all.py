import pytest
import os
import subprocess
import re
from pathlib import Path

import helpers


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
            (testcaserundir / laststatefile).symlink_to(os.path.relpath(base_directory / laststatefile, testcaserundir))
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
        helpers.assert_empty_stderr()
        helpers.assert_stdout_finished(message="GVEC SUCESSFULLY FINISHED!")
        helpers.assert_stdout_OpenMP_MPI()


@pytest.mark.post_stage
def test_post(
    runargs_prefix,
    binpath,
    testgroup,
    testcase,
    testcasepostdir,
    dryrun,
    annotations,
    artifact_pages_path,
):
    """
    Post processing  of statefile(s) from an example GVEC run.

    Takes the last two statefiles in `{rundir}/{testgroup}/{testcase}`  and runs the modified parameter with the post_gvec in
    `{postdir}/{testgroup}/{testcase}`.
    Note: the `testgroup` fixture is contained in `testcasepostdir`, but is given additionally to the test function
    for better readability and proper parameter ordering for the pytest nodeID.
    """
    args = runargs_prefix + [binpath / "gvec_post", "parameter.ini"]
    # find all statefiles in directory

    # run gvec
    with helpers.chdir(testcasepostdir):
        states = [sd for sd in os.listdir(".") if "State" in sd and sd.endswith(".dat")]
        for statefile in states[len(states) - 2 : len(states)]:
            args.append(statefile)  # add the last two files
        if dryrun:
            if len(states) == 0:
                args.append("STATEFILES???!")
            with open("dryrun-post.txt", "w") as file:
                file.write(f"DRYRUN: execute:\n {args} \n")
            return
        assert (
            len(states) > 0
        ), f"no statefile for post found in directory {testcasepostdir}"
        # run gvec_post
        with open("stdout.txt", "w") as stdout:
            stdout.write(f"RUNNING: \n {args} \n")
        with open("stdout.txt", "a") as stdout, open("stderr.txt", "w") as stderr:
            subprocess.run(args, text=True, stdout=stdout, stderr=stderr)
        # add link to artifact (CI)
        for filename in ["stdout", "stderr"]:
            if pages_postdir := os.environ.get("CASENAME"):
                pages_postdir = f"CIpost_{pages_postdir}"
            else:
                pages_postdir = "."
            annotations["gvec-output"].append(
                dict(
                    external_link=dict(
                        label=f"{testgroup}/{testcase}/{filename}",
                        url=f"{artifact_pages_path}/{pages_postdir}/{testgroup}/{testcase}/{filename}.txt",
                    )
                )
            )
        # check if GVEC was successful
        helpers.assert_empty_stderr()
        helpers.assert_stdout_finished(message="GVEC POST FINISHED !")
        helpers.assert_stdout_OpenMP_MPI()


@pytest.mark.converter_stage
@pytest.mark.parametrize(
    "which_conv", ["to_gene", "to_jorek", "to_castor3d", "to_hopr"]
)
def test_converter(
    runargs_prefix,
    binpath,
    testgroup,
    testcase,
    which_conv,
    testcaseconvdir,
    dryrun,
    annotations,
    artifact_pages_path,
):
    """
    Post processing  of statefile(s) from an example GVEC run, using the compiled converters.

    Takes the last statefile in `{rundir}/{testgroup}/{testcase}`  and runs the modified parameter with the post_gvec in
    `{postdir}/{which_conv}/{testgroup}/{testcase}`.
    Note: the `testgroup` fixture is contained in `testcasepostdir`, but is given additionally to the test function
    for better readability and proper parameter ordering for the pytest nodeID.
    """
    conv_def = {
        "to_gene": dict(exec="test_gvec_to_gene", msg="TEST GVEC TO GENE"),
        "to_hopr": dict(exec="test_gvec_to_hopr", msg="TEST GVEC TO HOPR"),
        "to_castor3d": dict(
            exec="convert_gvec_to_castor3d",
            msg="CONVERT GVEC TO CASTOR3D",
            args=[
                ["--rpoints=7", "--polpoints=12", "--torpoints=8", "--sflcoord=0"],
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
            ],
            fixedargs=[
                ["gvec2castor3d_sfl0.dat"],
                ["gvec2castor3d_sfl1.dat"],
                ["gvec2castor3d_sfl2.dat"],
                ["gvec2castor3d_sfl3.dat"],
            ],
        ),  # same length of args & fixedargs, give the number of runs
        "to_jorek": dict(
            exec="convert_gvec_to_jorek",
            msg="CONVERT GVEC TO JOREK",
            args=[["--rpoints=8", "--npfactor=1", "--polpoints=12"]],
            fixedargs=[["gvec2jorek_out.dat"]],
        ),
    }
    conv = conv_def[which_conv]
    if not ((not dryrun) and (binpath / conv["exec"]).exists()):
        pytest.skip(f"Executable {conv['exec']} not found in binary folder!")
        return
    # multiple runs with different arguments:
    if "args" in conv.keys():
        nruns = len(conv["args"])
        if "fixedargs" in conv.keys():
            assert len(conv["fixedargs"]) == len(conv["args"])
    else:
        nruns = 1
    # run converter
    with helpers.chdir(testcaseconvdir):
        for irun in range(0, nruns):
            args = runargs_prefix + [binpath / conv["exec"]]
            if "args" in conv.keys():
                args += conv["args"][irun]
            # find all statefiles in directory
            states = [
                sd for sd in os.listdir(".") if "State" in sd and sd.endswith(".dat")
            ]
            if dryrun:
                if len(states) == 0:
                    args.append("STATEFILES???!")
                    if "fixedargs" in conv.keys():
                        args.append(conv["fixedargs"][irun])
                with open(f"dryrun-post-converter{irun}.txt", "w") as file:
                    file.write(f"DRYRUN: execute:\n {args} \n")
                return
            assert (
                len(states) > 0
            ), f"no statefile for post-converter found in directory {testcaseconvdir}"

            args.append(states[-1])  # add the last state file
            if "fixedargs" in conv.keys():
                args += conv["fixedargs"][irun]
            # run gvec_post
            with open(f"stdout{irun}.txt", "w") as stdout:
                stdout.write(f"RUNNING: \n {args} \n")
            with (
                open(f"stdout{irun}.txt", "a") as stdout,
                open(f"stderr{irun}.txt", "w") as stderr,
            ):
                subprocess.run(args, text=True, stdout=stdout, stderr=stderr)
            # add link to artifact (CI)
            for filename in [f"stdout{irun}", f"stderr{irun}"]:
                if pages_convdir := os.environ.get("CASENAME"):
                    pages_convdir = f"CIconv_{pages_convdir}"
                else:
                    pages_convdir = "."
                annotations["gvec-output"].append(
                    dict(
                        external_link=dict(
                            label=f"{which_conv}/{testgroup}/{testcase}/{filename}",
                            url=f"{artifact_pages_path}/{pages_convdir}/{which_conv}/{testgroup}/{testcase}/{filename}.txt",
                        )
                    )
                )
            # check if GVEC was successful
            helpers.assert_empty_stderr(f"stderr{irun}.txt")
            helpers.assert_stdout_finished(
                f"stdout{irun}.txt", message=conv["msg"] + " FINISHED!"
            )
            # helpers.assert_stdout_OpenMP_MPI()


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
        if result in ["ignored", "symlink"]:
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
