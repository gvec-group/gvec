# Getting Started

## Virtual Environments

We strongly recommend you to use python's virtual environments to manage your python packages.
Virtual environments allow you to isolate the python packages you need for each project and prevent conflicts between the dependencies of different projects.
One way of creating and virtual environments using `venv` is given below:

1. make sure that you have `venv` installed, e.g. on Ubuntu the `python3-venv` apt-package
2. create a new virtual environment in your current directory with
   ```bash
   python3 -m venv .venv
   ```
   This will now create a `.venv` directory which can contain python packages and executables
3. activate the environment
   ```bash
   source .venv/bin/activate
   ```
   This will now cause `python` and `pip` to use this environment, as well as make any executables in that environment available in your `PATH`.

## Installing GVEC

In the best case, GVEC can be installed with its python bindings using:
```bash
pip install gvec
```
For system requirements and more information see the [installation](/user/install) page.

To confirm GVEC is installed and working, call `pygvec --version` in the terminal, or `import gvec` in a python session.

## On the tutorials

Each subsection is a Jupyter notebook, displayed after its execution.
You can follow the link in the right sidebar to download or browse the source code.

We recommend going through the notebooks in the order of this section.

You can also simply [download all notebooks](<path:/static/gvec_tutorials.zip>) as a zip file.

The available tutorials are:
   * [](/tutorials/notebooks/010_tokamak)
   * [](/tutorials/notebooks/020_stellarator)
   * [](/tutorials/notebooks/021_boozer)
   * [](/tutorials/notebooks/030_stages)
   * [](/tutorials/notebooks/040_current)
   * [](/tutorials/notebooks/050_pygvec)
   * [](/tutorials/notebooks/060_coils)
   * [](/tutorials/notebooks/061_coils)

### Running Jupyter notebooks

The GVEC tutorials are mostly organized in *Jupyter notebooks* (also called *ipython notebooks*).
You can however also follow along with a normal python script.

To use the notebooks you need to have the following packages installed:
```bash
pip install jupyter ipykernel ipympl ipywidgets matplotlib
```

We recommend two options for running Jupyter notebooks:
1. Visual Studio Code with the `jupyter` extension
2. The jupyter notebook browser interface

   Here need to make the virtual environment visible (`.venv` must the activated!) with
   ```
   python3 -m ipykernel install --user --name .venv
   ```
   and then start
   ```bash
   jupyter notebook &
   ```

When you **open a notebook**, you should be prompted which kernel to use.
Select the environment where you installed `gvec`, e.g. `.venv`.
