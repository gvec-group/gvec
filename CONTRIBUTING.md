# Contributing to GVEC

## Generating Documentation

* The code documentation is auto-generated with [FORD](https://forddocs.readthedocs.io/en/latest/), which can be installed with `pip install ford`
* Static documentation (User / Developer Guides) are written in markdown and defined in `docs/static`
* To generate the documentation run `ford docs/ford.md` from within the main repository. This will generate a `docs/_build` folder containing all the documentation in HTML & CSS.
* A job within the CI script generates the documentation for each `develop` pipeline and deploys it to GitLab pages.