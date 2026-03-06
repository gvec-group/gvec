---
file_format: mystnb
kernelspec:
  name: python3
---
# pygvec CLI

Users can interface with `pygvec` directly from the command line rather than by opening a python REPL or Jupyter window:

```{code-cell} ipython3
:tags: [remove_input]
!pygvec --help
```

To get help for a specific subcommand, use `pygvec <subcommand> --help`.

Most of the subcommands also support different verbosity levels.
Per default, warnings and errors are printed, but are suppressed with `-q`/`--quiet`.
The first verbosity level (`-v`/`--verbose`) shows additional *INFO* messages.
The second verbosity level (`-vv`) also shows *DEBUG* messages, and `pygvec run` supports a third verbosity level (`-vvv`) that includes the full GVEC diagnostic output which is otherwise redirected to a logfile.

The `pygvec` interfaces for other codes or libraries are:
- [quasr](./quasr.md)
- [cas3d](./cas3d.md)
- [gist](./gist.md)

There are also tools for visualisation with:
- `pygvec visu`: for generating some default plots. Users should see the [plotting](/tutorials/notebooks/051_plotting) tutorial for a comprehensive guide for using the plotting routines.
- `pygvec to-paraview`: for exporting a GVEC run to a paraview file for visualisation.
