# Developer's Notes

## Building documentation
The documentation of GVEC is split into three parts:
1) The mathematical details, found in the *GVEC Prototype* pdf
2) The user & developer documentation written in *restructured text* and *markdown* and compiled with *sphinx*
3) Auto-generated fortran code documentation using *ford*

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

The required python packages to build the documentation are given in `docs/requirements.txt` and can be installed with:
```bash
pip install -r docs/requirements.txt
```
In addition the `graphs` feature of *ford* requires an installation of *graphviz/dot* and *urw-fonts*.