# Developer's Notes

```{toctree}
:titlesonly:
testing.md
pipeline.md
```

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

## Object-Oriented Programming in FORTRAN

Here is a recommendation for a tutorial on how to program in an object-oriented way
with [polymorphism in fortran](https://gist.github.com/n-s-k/522f2669979ed6d0582b8e80cf6c95fd).

## Testing

*   GVEC has been equipped with automatic CI testing on the `gitlab.mpcdf.mpg.de` server, using shared MPCDF gitlab runners to execute the tests. 
    *   More details on the CI setup are found at <project:pipeline.md>.
    *   The CI manages different builds of the code, then calls pytest for running them and checking the results (requires `python >3.10` to be installed!).
*   The `pytest` feature also **allows to locally run** the same tests. More details and examples on running the tests with pytest are found at <project:testing.md>.
*   A predefined set of tests can be executed using `ctest`, after the [cmake install process](project:INSTALL.md). Simply change to the build directory, and execute:
    ```bash
    ctest -T test --output-on-failure -R
    ```