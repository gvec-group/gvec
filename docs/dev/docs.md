# Documentation

## Writing documentation
* add *markdown* or *restructured text* files in the `docs` directory and link them in the TOC of `docs/index.rst`
* files & directories in `docs/`:
    * `index.rst` is the landing page and contains the table of contents
    * `conf.py` contains the *sphinx* configuration as a python script
        * we use the [pydata-sphinx-theme](https://pydata-sphinx-theme.readthedocs.io)
        * and the [myst-parser](https://myst-parser.readthedocs.io) for extended *markdown*
    * `Makefile` contains the logic to build the *sphinx* documentation
    * `requirements.txt` contain the python packages required to build the documentation (see [](#building-documentation))
    * `templates/` contains html templates that can be used to style the documentation
    * `static/` contains content that should be copied to the build directory
    * `ford/` contains the [ford](https://forddocs.readthedocs.io) configuration (`ford.md`) and auxiliary files
    * `prebuild/` contains auto-generated content (e.g by *ford*) to be included in the build directory
    * `build/` is the default directory where the documentation output is saved to, e.g. within `build/html/`
* build the documentation locally and run a manually triggered *scheduled pipeline*

## Building documentation
The documentation of GVEC is split into three parts:
1) Mathematical details, found in [theory and implementation details](https://gitlab.mpcdf.mpg.de/gvec-group/GVEC_doc/blob/master/GVEC_prototype/GVEC_prototype.pdf) 
2) [User and developer documentation](https://gvec-group.pages.mpcdf.de/gvec) written in *restructured text* and *markdown* and compiled with [sphinx](https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html)
3) Auto-generated [fortran code documentation](https://gvec-group.pages.mpcdf.de/gvec/_static/ford/index.html) built with [ford](https://forddocs.readthedocs.io)

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

The required python packages to build the documentation are given in `docs/requirements.txt` and can be installed with:
```bash
pip install -r docs/requirements.txt
```
In addition the `graphs` feature of *ford* requires an installation of *graphviz/dot* and *urw-fonts*.
