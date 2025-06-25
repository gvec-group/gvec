# Developer Guide

These pages contain development guidelines and useful resources for developing GVEC.

## Contact

GVEC is mainly being developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
led by Prof. Eric Sonnendruecker at the Max Planck Institute for Plasma Physics
in Garching, Germany. Outside contributions are of course very welcome!

<!-- Other Topics -->

## Development Workflow

* main repository on [MPCDF-GitLab](https://gitlab.mpcdf.mpg.de/gvec-group/gvec): Issues, Merge Requests, CI-Pipelines
    * requires a MPCDF account, contact the maintainers to obtain a guest account
* mirror repository on [GitHub](https://github.com/gvec-group/gvec) for visibility and public contributions
* prefer merging over rebasing
* automatic testing of all pushes to GitLab, **add tests for new features**
* use feature branches, merge to `develop` early and often (*at least in theory*)
    * use GitLab *merge requests* to document the changes and code review
* `main` points to the latest release / tag
    * releases (with corresponding tags) are created within GitLab
    * associate milestones with the releases to document progress
    * tags/releases (mostly) follow [semantic versioning](https://semver.org/)
* use pre-commit hooks (python formatting, notebook cleaning, etc.)
    * `pip install pre-commit` & `pre-commit install`

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
Contributors <CONTRIBUTORS>
```

```{toctree}
:caption: API
gvec.core.state <api/core-state>
gvec.core.run <api/core-run>
gvec.core.compute <api/core-compute>
gvec.quantities <api/quantities>
gvec.fourier <api/fourier>
gvec.surface <api/surface>
gvec.util <api/util>
gvec.vtk <api/vtk>
gvec.scripts.main <api/scripts-main>
gvec.scripts.run <api/scripts-run>
gvec.scripts.cas3d <api/scripts-cas3d>
gvec.scripts.quasr <api/scripts-quasr>
gvec.lib <api/lib>
```
