# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
# with inspiration from https://gitlab.mpcdf.mpg.de/struphy/struphy/-/blob/devel/doc/conf.py

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import subprocess
import logging


# -- Project information -----------------------------------------------------

project = 'GVEC'
copyright = '2017-2024 (c) Florian Hindenlang | Max Planck Institute for Plasma Physics'
author = 'Florian Hindenlang et al. | Max Planck Institute for Plasma Physics'

try:
    p = subprocess.run(["git", "describe", "--tags", "--dirty", "--always"], capture_output=True)
    version = p.stdout.decode().strip()
except Exception as e:
    logging.error(f"Could not get git version: {e}")
    version = "unknown"
release = version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["ford/ford.md"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

html_theme_options = {
    # from struphy:
    "sidebarwidth": 270,
    "show_nav_level": 3,
    "show_toc_level": 2,
    "navigation_depth": 4,
    "header_links_before_dropdown": 7,
    "primary_sidebar_end": ["sidebar-ethical-ads"],

    "footer_start": ["copyright", "version", "last-updated"],
    "footer_end": ["sphinx-version", "theme-version"],

    "gitlab_url": "https://gitlab.mpcdf.mpg.de/gvec-group/gvec",
    "external_links": [
        {"name": "Fortran Code Documentation", "url": "_static/ford/index.html"},
    ]
}

html_sidebars = {
   '**': ['globaltoc.html', 'relations.html', 'searchbox.html'],
}

html_title = "GVEC Documentation"
html_last_updated_fmt = "%Y-%m-%d"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['static', 'prebuild']
# static-build contains pages generated outside of sphinx - e.g. ford
# static is intended for truly static pages

# --- markdown parsing with myst --- #
myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

myst_dmath_allow_labels=True
