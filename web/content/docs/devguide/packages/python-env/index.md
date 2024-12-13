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

- [TESPy]({{< ref "3d_3bhes_array.md#tespy" >}}) for simulating thermal engineering plants in a benchmark
- [pvpython](https://www.paraview.org/paraview-docs/latest/python/) for pre- and post-processing
- [nbconvert](https://nbconvert.readthedocs.io/en/latest/) for testing Jupyter Notebooks
- ...

Python packages are usually installed via `pip` inside an isolated environment (a virtual environment).

## Pip

When configuring OGS with `OGS_USE_PIP=ON` Python creates a new virtual environment in the `.venv`-directory inside your build directory. It will also install required Python packages into this environment. You can see the current environment definition in the file `requirements.txt` inside your build-directory.

Make sure, that you did

```bash
sudo apt-get install python3 python3-pip
sudo apt-get install python3-venv
```

then you can use the `-DOGS_USE_PIP=ON` option:

```bash
cmake ../ogs -DOGS_USE_PIP=ON -DCMAKE_BUILD_TYPE="Release" -G Ninja
```

When configuring OGS with `OGS_USE_PIP=ON`, Python creates a new virtual
environment in the .venv-directory inside your build directory.
It will also install required Python packages into this environment.
For example, ogstools or Jupyter will this way be available ready-made in a
consistent configuration.
You can see the current environment definition with all the included packages in
the file requirements.txt inside your build-directory.

When you want to use the python packages in this virtual environment, you need
to activate it by

```bash
source <your-build-dir>/.envrc
```

where `.envrc` is a little script including `source .venv/bin/activate` as well as some path settings.

To manually add Python packages run the following inside your build-directory:

```bash
.venv/bin/pip install python-package-name
```

To activate the environment run `source .envrc` inside your build directory. If you have the [`direnv`](https://direnv.net)-tool installed and setup the virtual environment will be activated automatically upon changing into the build directory.

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
