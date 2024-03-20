import pytest
import sys, os
import subprocess
import re
from pathlib import Path
import shutil

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


@pytest.fixture(scope="session")
def testcaserundir(rundir: Path, testgroup: str, testcase: str):
    """
    Generate the run directory at `{rundir}/{testgroup}/{testcase}` based on `{testgroup}/{testcase}`
    """
    # assert that `{rundir}` and `{rundir}/{testgroup}` exist
    if not rundir.exists():
        rundir.mkdir()
    if not (rundir / "data").exists():
        (rundir / "data").symlink_to(Path(__file__).parent / "data")
    if not (rundir / testgroup).exists():
        (rundir / testgroup).mkdir()
    # create the testcase directory
    sourcedir = Path(__file__).parent / "examples" / testcase
    targetdir = rundir / testgroup / testcase
    if targetdir.exists():
        shutil.rmtree(targetdir)
    shutil.copytree(sourcedir, targetdir, symlinks=True)
    # special treatment
    match testgroup:
        case "shortrun":
            helpers.adapt_parameter_file(
                sourcedir / "parameter.ini",
                targetdir / "parameter.ini",
                testlevel="-1",
                logIter="1",
                MaxIter="1",
                outputIter="1",
            )
        case "debugrun":
            helpers.adapt_parameter_file(
                sourcedir / "parameter.ini",
                targetdir / "parameter.ini",
                testlevel="2",
                logIter="1",
                MaxIter='1',
                outputIter="1",
            )
    return targetdir


@pytest.fixture(scope="session")
def testcasepostdir(postdir: Path, rundir: Path, testgroup: str, testcase: str):
    """
    Generate the post directory at `{postdir}/{testgroup}/{testcase}` based on `{rundir}/{testgroup}/{testcase}`
    """
    # assert that `{postdir}` and `{postdir}/{testgroup}` exist
    if not postdir.exists():
        postdir.mkdir()
    if not (postdir / "data").exists():
        (postdir / "data").symlink_to(Path(__file__).parent / "data")
    if not (postdir / testgroup).exists():
        (postdir / testgroup).mkdir()
    # create the testcase directory
    sourcedir = Path(__file__).parent / "examples" / testcase
    sourcerundir = rundir / testgroup / testcase
    targetdir = postdir / testgroup / testcase
    if targetdir.exists():
        shutil.rmtree(targetdir)
    # copy input files from examples/testcase
    shutil.copytree(sourcedir, targetdir, symlinks=True)
    states = [
        sd for sd in os.listdir(sourcerundir) if "State" in sd and sd.endswith(".dat")
    ]
    # link to statefiles from run_stage
    for statefile in states:
        (targetdir / statefile).symlink_to(sourcerundir / statefile)
    # overwrite parameter file with the rundir version and modify it
    helpers.adapt_parameter_file(
        sourcerundir / "parameter.ini",
        targetdir / "parameter.ini",
        visu1D="!0",
        visu2D="!0",
        visu3D="!0", 
        SFLout="!0",  # only uncomment visualization flags
    )
    return targetdir


# === TESTS === #


@pytest.mark.run_stage
def test_run(binpath, testgroup, testcaserundir, testcase, dryrun, annotations, artifact_pages_path):
    """
    Test end2end GVEC runs with `{testgroup}/{testcase}/parameter.ini`

    Note: the `testgroup` fixture is contained in `testcaserundir`, but is given additionally to the test function
    for better readability and proper parameter ordering for the pytest nodeID.
    """
    args = [binpath / "gvec", "parameter.ini"]
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
            (testcaserundir / laststatefile).symlink_to(base_directory / laststatefile)
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
            annotations["gvec-output"].append(dict(external_link=dict(
                label=f"{testgroup}/{testcase}/{filename}", 
                url=f"{artifact_pages_path}/{pages_rundir}/{testgroup}/{testcase}/{filename}.txt")))
        # check if GVEC was successful
        helpers.assert_empty_stderr()
        helpers.assert_stdout_finished(message="GVEC SUCESSFULLY FINISHED!")


@pytest.mark.post_stage
def test_post(binpath, testgroup, testcase, testcasepostdir, dryrun):
    """
    Post processing  of statefile(s) from an example GVEC run.

    Takes the last two statefiles in `{rundir}/{testgroup}/{testcase}`  and runs the modified parameter with the post_gvec in
    `{postdir}/{testgroup}/{testcase}`.
    Note: the `testgroup` fixture is contained in `testcasepostdir`, but is given additionally to the test function
    for better readability and proper parameter ordering for the pytest nodeID.
    """
    args = [binpath / "gvec_post", "parameter.ini"]
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
        # check if GVEC was successful
        helpers.assert_empty_stderr()
        helpers.assert_stdout_finished(message="GVEC POST FINISHED !")


@pytest.mark.regression_stage
def test_regression(testgroup, testcase, rundir, refdir, dryrun, logger, reg_rtol, reg_atol, extra_ignore_patterns):
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
        pytest.skip(f"Testcase does not exist")
    if not testcaserefdir.exists():
        pytest.skip(f"Reference does not exist")
    # compare the list of files in the two directories
    runfiles, reffiles = (
        set([file for file in os.listdir(directory) if "dryrun" not in file and not file.endswith("~")])
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
    num_diff_files = 0
    num_diff_lines = 0
    num_warnings = 0
    for filename in runfiles & reffiles:
        # skip symbolic links (mostly unreachable input)
        if (testcaserundir / filename).is_symlink():
            results[filename] = "symlink"
            continue
        # statefiles
        elif re.match(r"\w+State\w+\.dat", filename):
            num = helpers.check_diff_files(
                testcaserundir / filename,
                testcaserefdir / filename,
                ignore_regexs=[r".*/.*"] + extra_ignore_patterns, # ignore lines with a path
                atol=reg_atol,
                rtol=reg_rtol,
            )
        elif re.match(r"log\w*\.csv", filename):
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
                ignore_regexs=[r".*GIT_.*",r".*CMAKE.*",r".*sec.*", r".*date.*", r".*PosixPath.*", r"^[\s=]*$"] + extra_ignore_patterns,
                warn_regexs=["Number of OpenMP threads"],
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
    if num_warnings > 0:
        logger.warning(f"Found {num_warnings} warnings!")
        pytest.raised_warnings = True
    if num_diff_files > 0 or runfiles != reffiles:
        msg = f"Found {num_diff_files} different files with {num_diff_lines} different lines, " \
            f"{len(runfiles - reffiles)} additional files and {len(reffiles - runfiles)} missing files."
        logger.info(f"{' SUMMARY ':=^80}")
        logger.error(msg)
        red, reset = "\x1b[31;20m", "\x1b[0m"
        for filename, result in sorted(results.items(), key=lambda x: (x[1], x[0])):
            if result in ["success", "ignored", "symlink"]:
                logger.debug(f"... {result} : {filename}")
            else:
                logger.error(f"... {red}{result}{reset} : {filename}")
        for filename in runfiles - reffiles:
            logger.error(f"... {red}extra{reset}   : {filename}")
        for filename in reffiles - runfiles:
            logger.error(f"... {red}missing{reset} : {filename}")
        raise AssertionError(msg)
    
