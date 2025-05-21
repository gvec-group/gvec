# Documentation

The documentation of GVEC is split into the following parts:
1) [User and developer documentation](/index) written in *restructured text* and *markdown* and compiled with [sphinx](https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html)
2) Auto-generated [fortran code documentation](../ford/index.html){.external} built with [ford](https://forddocs.readthedocs.io)

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
    * `generators/` contains additional material that generates figures. And there is a script for generating the parameter list (`generate_parameter_list.py`), which parser the `parameters.yaml` file. This is called in `conf.py` to generate markdown files `parameters-*.md`. To add a new parameter to the yaml file, the `param_dict_to_doc.ipynb` has to be edited and run.
    * `ford/` contains the [ford](https://forddocs.readthedocs.io) configuration (`ford.md`) and auxiliary files
    * `ford/static` contains the static pages processed by *ford*, currently only a redirect to the main documentation is used.
    * `extra/` contains auto-generated content (e.g by *ford*) to be included in the build directory
    * `build/` is the default directory where the documentation output is saved to, e.g. within `build/html/`
* you can build the documentation locally, run a manually triggered *scheduled pipeline* or manually run a *publish/pages* job

## Building documentation

The documentation is built and deployed on *Read the Docs* with the url [https://gvec.readthedocs.io](https://gvec.readthedocs.io).
*Read the Docs* allows switching between different versions.
Which version are built is configured in the [Read the Docs Settings](https://app.readthedocs.org/projects/gvec).
Supported versions users can switch between are configured in `docs/static/version-switcher.json`.
*Read the Docs* also provides a small menu in the lower right corner that can be used to switch versions, e.g. to test a branch with a new version.

The current strategy is to provide the documentation for the latest release (the `main` branch) under [https://gvec.readthedocs.io/main](https://gvec.readthedocs.io/main)
and the contents of the `docs` branch under the default [https://gvec.readthedocs.io/latest](https://gvec.readthedocs.io/latest).
The `docs` branch should be kept up to date with `develop` and merged frequently, but allows us to commit directly to it, updating the documentation with minimal overhead.

[FORD](https://forddocs.readthedocs.io/en/latest/) is configured in the `docs/ford/ford.md` file and can manually be triggered with:
```bash
ford docs/ford/ford.md
```
This will generate files in `docs/extra/ford`.

The sphinx documentation is configured in `docs/conf.py` and `docs/Makefile` and build with
```bash
cd docs
make cleanall
make html
```
generating documentation in `docs/build`. It will copy files from `docs/static` and `docs/extra`.
If your webbrowser cannot render the html and css files in `docs/build/html/` properly, you can start a local webserver with `python -m http.server`

### Internals

The `FORD_PREFIX` environment variable is used to set a path prefix to the absolute path used to link to the ford documentation in the header. By default this is `/ford/index.html`, but on GitLab we need to link to `/gvec/ford/index.html`.
