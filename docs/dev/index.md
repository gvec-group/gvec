# Developer Guide

These pages contain development guidelines and useful resources for developing GVEC.

## Contact

GVEC is mainly being developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
led by Prof. Eric Sonnendruecker at the Max Planck Institute for Plasma Physics
in Garching, Germany. Outside contributions are of course very welcome!

If you have questions you can contact the maintainers ([Florian](https://www.ipp.mpg.de/person/139885/5497858)) directly or join the [matrix-chat channel](https://matrix.to/#/#gvec:mpg.de).

<!-- Other Topics -->

## Development Workflow

GVEC is hosted on two *git* repositories:
1) the main development repository on [MPCDF-GitLab](https://gitlab.mpcdf.mpg.de/gvec-group/gvec): Issues, Merge Requests, CI-Pipelines
    * Benefit: in-house support, unlimited CI/CD budget
    * Drawback: requires an MPCDF account for contributions → contact the maintainers to obtain a guest account
2) a public repository on [GitHub](https://github.com/gvec-group/gvec) for visibility and public contributions

### New contributors

**Contributions are always welcome!**

If you already have a GitHub account, just [open an issue](https://github.com/gvec-group/gvec/issues/new) on GitHub if you want to report a bug or have a question.

If you want to contribute, it is easiest to submit a Pull Request on GitHub for your proposed changes. The Maintainers will then mirror your branch back to GitLab to run all the testing pipelines. For larger contributions, you can contact the maintainers to obtain a guest account for the MPCDF-GitLab.

You can also join the [matrix-chat channel](https://matrix.to/#/#gvec:mpg.de) for questions and discussions.

### Remarks / Guidelines

* `develop` is the main development branch
    * use feature branches, merge to `develop` early and often
    * use GitLab *merge requests* to document the changes and code review
    * update with `develop` before merging → resolve conflicts in the feature branch
* prefer merging over rebasing
* automatic testing in GitLab: **add tests for new features**
* `main` points to the latest release / tag
    * releases (with corresponding tags) are created within GitLab
    * associate milestones with the releases to document progress
    * tags/releases (mostly) follow [semantic versioning](https://semver.org/)
* fixes should be added to `develop` / a feature branch and can be cherry-picked back to older releases if necessary
    * a release branch (e.g. `v1.3.x`) can be setup to track the history of a release which has branched away from `develop`
* `docs` can be used to quickly fix the documentation (allows directly committing to it)
    * On readthedocs `latest` → `docs` branch, the other branches / tags are used by name
    * configure readthedocs here: https://app.readthedocs.org/projects/gvec/
    * in repo: `.readthedocs.yaml` & `docs/static/version-switcher.json`
* use pre-commit hooks (python formatting, notebook cleaning, etc.)
    * `pip install pre-commit` & `pre-commit install`
    * use *ruff* for python formatting & linting
* split your work into small commits
* don't commit large binary data (use pre-commit!)

## Repository structure

* `src/` - the main fortran sources
* `pyproject.toml` - configuration for the `gvec` python package
* `python/gvec/` - the `gvec` python package (pyGVEC) with bindings to fortran
* `python/examples/` - example notebooks and configurations
* `python/kind_map.py` & `python/class_names.py` - auxiliary files for the python bindings with *f90wrap*
* `CMakeLists.txt`, `CMakePresets.json`, `cmake/` - configuration of *CMake*
* `CI_setup/` - scripts to load modules for different clusters & CI runners
* `test-CI/` - testcases and test logic using `pytest`
* `.gitlab-ci.yml` & `CI_templates` - configuration of the GitLab CI Pipelines (see <dev/pipeline>)
* `docs/` & `.readthedocs.yaml` - configuration and static content for the documentation, built with *sphinx* and *ford*
* `.gitignore` - file patterns to be ignored by *git*
* `.mailmap` - cleaning git authors for `git blame`
* `template/` - a structural template for fortran sources
* `tools/`

## Object-Oriented Programming in FORTRAN

Here is a recommendation for a tutorial on how to program in an object-oriented way
with [polymorphism in fortran](https://gist.github.com/n-s-k/522f2669979ed6d0582b8e80cf6c95fd).

## Useful VSCode extensions

* Modern Fortran
* CMake Tools
* Git Graph
* GitLab Workflow
* GitLens (Premium/Students)
* GitHub Copilot (AI, Premium/Students)
* Codeium (AI)
* Jupyter
* MyST-Markdown
* Python
* Pylance
* Ruff (Python Linter & Formatter)
* Todo Tree
* Vim
* YAML
* netCDF Preview

## Contents

<!-- TOC -->

```{toctree}
:caption: Developer Guide

testing
pipeline
docs
python
fortran
derivations
Contributors <CONTRIBUTORS>
```

```{toctree}
:caption: API
gvec.core.state </api/core-state>
gvec.core.run </api/core-run>
gvec.core.compute </api/core-compute>
gvec.quantities </api/quantities>
gvec.fourier </api/fourier>
gvec.surface </api/surface>
gvec.util </api/util>
gvec.vtk </api/vtk>
gvec.gframe </api/gframe>
gvec.scripts.main </api/scripts-main>
gvec.scripts.run </api/scripts-run>
gvec.scripts.cas3d </api/scripts-cas3d>
gvec.scripts.quasr </api/scripts-quasr>
gvec.lib </api/lib>
```
