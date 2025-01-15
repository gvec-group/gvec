# Developer Guide

These pages contain development guidelines and useful resources for developing GVEC.

## Contact

GVEC is mainly being developed in the department of **Numerical Methods in Plasma Physics (NMPP)**
led by Prof. Eric Sonnendruecker at the Max Planck Institute for Plasma Physics
in Garching, Germany. Outside contributions are of course very welcome!

<!-- TOC -->

```{toctree}
:caption: Developer Guide

testing
pipeline
docs
Contributors <CONTRIBUTORS>
```

<!-- Other Topics -->

## Development Workflow

* prefer merging over rebasing
* use `git-lfs` for large files
* automatic testing of all pushes to GitLab, **add tests for new features**
* use feature branches, merge to `develop` early and often (*at least in theory*)
    * use GitLab *merge requests* to document the changes and code review
* `main` points to the latest release / tag
    * releases (with corresponding tags) are created within GitLab
    * associate milestones with the releases to document progress

## Repository structure

* `src/` - the main fortran sources
* `ini/` - example parameter files for various configurations
* `CMakeLists.txt` & `CMakePresets.json` - configuration of CMake
* `cmake/` - additional CMake configuration for converters
* `CI_setup/` - scripts to load modules for different clusters & CI runners
* `test-CI/` - testcases and test logic using `pytest`
* `.gitlab-ci.yml` & `CI_templates` - configuration of the GitLab CI Pipelines (see <dev/pipeline>)
* `docs/` - configuration and static content for the documentation
* `.gitignore` - file patterns to be ignored by `git`
* `.gitattributes` - configuration of `git-lfs`
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
* ~~Black Formatter~~
* Ruff (Python Linter & Formatter)
* Todo Tree
* Vim
* YAML
* netCDF Preview
