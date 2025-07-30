# Introduction

The tutorials in this section are organized in ipython notebooks.
Each subsection is generated from one notebook, displayed after its execution. You can follow the link in the right sidebar to download or browse the source code.

## Requirements

The requirements to execute the ipython notebooks locally are:
1. have [gvec and its python bindings installed](install.md) in a virtual environment (`.venv`)
1. Additionally install the following packages in the same environment:

   `matplotlib ipywidgets ipykernel ipympl jupyter`
1. you need to be able to run an `ipython` notebook (`.ipynb` file).
   We tested two options:
   1. VScode with `jupyter` extension installed
   2. jupyter notebook via a browser.
      Here need to make the virtual environment visible (`.venv` must the activated!) with
      ```
      python3 -m ipykernel install --user --name .venv
      ```
      and then start
      ```bash
      jupyter notebook &
      ```

When you open a notebook, you should be prompted which kernel to use.
**Be sure** to select the virtual environment `.venv` as a kernel!
