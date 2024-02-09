# Automatic test system

## Usage
* Call `pytest` or `python -m pytest` in the `test-CI` directory to run automatic tests.
* With the `-m MARKERS` option tests can be selected or deselected based on their markers.
    * Currently the three custom markers are: `example`, `shortrun` and `regression`.
    * The default is `-m "not regression"`.
    * Example: `-m example` selects only tests marked with example, `-m "not example"` all but those tests
* The relevant paths are supplied to `pytest` with `--builddir`, `--rundir` and `--refdir`.

### Examples
* `python -m pytest -v --rundir RUNDIR --builddir BUILDDIR` to run all end2end tests using `BUILDDIR/bin/gvec` and store the results at `RUNDIR`
    * omitting `--builddir BUILDDIR` sets the default `BUILDDIR=build`
    * omitting `--rundir RUNDIR` sets the default `RUNDIR=test-CI/run`
* `python -m pytest -v -m shortrun` to only run shortrun tests
* `python -m pytest -v -m regression --rundir RUNDIR --refdir REFDIR` to run the regression tests by comparing RUNDIR to REFDIR

## Details
* Large input files should be stored in `test-CI/data` and the examples should contain a *relative* symbolic link
    * When filling a `RUNDIR`, a link `RUNDIR/data -> test-CI/data` is setup. In this way the relative link should still work.
* the test logic is contained in `test_examples.py:test_examples` but relies on other functions:
    * the test logic for a successful GVEC run is in `helpers.py:assert_empty_stderr` and `assert_finished_stdout`
    * the test logic for statefile comparison in `helpers.py:assert_equal_statefiles`