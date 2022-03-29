+++
date = "2020-10-05T15:16:13+01:00"
title = "Python environment"
author = "Lars Bilke"
weight = 1040

aliases = ["/docs/devguide/advanced/python-env"]

[menu]
  [menu.devguide]
    parent = "packages"
+++

In OGS we make use of Python packages at different stages, e.g.:

- [ogs-container-maker](https://gitlab.opengeosys.org/ogs/container-maker) when the CI prepares its environment
- [TESPy]({{< ref "3d_3bhes_array.md#tespy" >}}) for simulating thermal engineering plants in a benchmark
- [pvpython](https://kitware.github.io/paraview-docs/latest/python/) for pre- and post-processing
- [nbconvert](https://nbconvert.readthedocs.io/en/latest/) for testing Jupyter Notebooks
- ...

Python packages are usually installed via `pip` inside an isolated environment (a virtual environment).

## Pip

When configuring OGS with `OGS_USE_PIP=ON` Python creates a new virtual environment in the `.venv`-directory inside your build directory. It will also install required Python packages into this environment. You can see the current environment definition in the file `requirements.txt` inside your build-directory.

To manually add Python packages run the following inside your build-directory:

```bash
.venv/bin/pip install python-package-name
```

To activate the environment run `source .venv/bin/activate` inside your build directory.

### Pip & Benchmarks

You can use the argument `PYTHON_PACKAGES` on `AddTest()` to specify additional Python package dependencies.

The following example would install the latest version of `numpy` and `pandas` version 0.1.2:

```cmake
AddTest(
    ...
    PYTHON_PACKAGES numpy pandas==0.1.2
    ...
)
```
