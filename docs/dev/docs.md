# Documentation

The documentation of GVEC is split into three parts:
1) Mathematical details, found in [theory and implementation details](https://gitlab.mpcdf.mpg.de/gvec-group/GVEC_doc/blob/master/GVEC_prototype/GVEC_prototype.pdf) 
2) [User and developer documentation](/index) written in *restructured text* and *markdown* and compiled with [sphinx](https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html)
3) Auto-generated [fortran code documentation](/_static/ford/index.html){.external} built with [ford](https://forddocs.readthedocs.io)

## Requirements
The required python packages to build the documentation are given in `docs/requirements.txt` and can be installed with:
```bash
pip install -r docs/requirements.txt
```
<!-- Apparently the consensus in the python community is to keep the development dependencies (e.g. for building the docs) in a seperate `requirements.txt` file and not in `pyproject.toml`. -->

In addition the `graphs` feature of *ford* requires an installation of *graphviz/dot* and *urw-fonts*.


## Writing User & Developer documentation
* static pages (guides, examples, etc.) are written in *restructured text* or *markdown*
    * we use the [myst-parser](https://myst-parser.readthedocs.io) for an extended markdown supporting most sphinx directives
* content is grouped into two subdirectories: `docs/user` for user documentation and `doc/dev` for developer documentation
    * each directory corresponds to a section in the documentation, i.e. different left sidebar for navigation
* the third section is the auto-generated fortran api using [ford](https://forddocs.readthedocs.io)
* add all new pages to the *toctree* in the respective `index.md`
* files & directories in `docs/`:
    * `index.rst` is the landing page and contains the main table of contents
    * `conf.py` contains the *sphinx* configuration as a python script
        * please add comments when you extend the configuration
    * `Makefile` contains the logic to build the *sphinx* documentation
    * `requirements.txt` contain the python packages required to build the documentation
    * `templates/` contains html templates that can be used to style the documentation
    * `static/` contains content that should be copied directly to the build directory
    * `ford/` contains the [ford](https://forddocs.readthedocs.io) configuration (`ford.md`) and auxiliary files
    * `ford/static` contains the static pages processed by *ford*, currently only a redirect to the main documentation is used.
    * `prebuild/` contains auto-generated content (e.g by *ford*) to be included in the build directory
    * `build/` is the default directory where the documentation output is saved to, e.g. within `build/html/`
* you can build the documentation locally, run a manually triggered *scheduled pipeline* or manually run a *publish/pages* job

## Building documentation

Parts 2 and 3 are compiled using the CI and deployed to GitLab Pages. In the future we might also deploy to *readthedocs*.

[FORD](https://forddocs.readthedocs.io/en/latest/) is configured in the `docs/ford/ford.md` file and can manually be triggered with:
```bash
ford docs/ford/ford.md
```
This will generate files in `docs/prebuild/ford`.

The sphinx documentation is configured in `docs/conf.py` and `docs/Makefile` and build with
```bash
cd docs
make html
```
generating documentation in `docs/build`. It will copy files from `docs/static` and `docs/prebuild`.
If your webbrowser cannot render the html and css files in `docs/build/html/` properly, you can start a local webserver with `python -m http.server`
