+++
date = "2020-10-05T15:16:13+01:00"
title = "Python environment"
author = "Lars Bilke"
weight = 1030

aliases = ["/docs/devguide/advanced/python-env"]

[menu]
  [menu.devguide]
    parent = "packages"
+++

In OGS we make use of Python packages at different stages, e.g.:

- [Conan]({{< ref "conan.md" >}}) at configure-time to install third-party dependencies
- [ogs-container-maker](https://gitlab.opengeosys.org/ogs/container-maker) when the CI prepares its environment
- [TESPy]({{< ref "3d_3bhes_array.md#tespy" >}}) for simulating thermal engineering plants in a benchmark
- [pvpython](https://kitware.github.io/paraview-docs/latest/python/) for pre- and post-processing
- ...

Python packages are usually installed via `pip` and `virtualenv` provides an isolated environment to install these packages to.

## Poetry

We make use of [poetry](https://python-poetry.org) which is a wrapper for `pip` and `virtualenv`. We recommend to install `poetry` with your system package manager or:

<div class='win'>

```ps
(Invoke-WebRequest -Uri https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py -UseBasicParsing).Content | python -
```

</div>

<div class='linux'>

```bash
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
```

</div>

<div class='mac'>

```bash
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
```

</div>

When configuring with CMake `poetry` creates a new virtual environment in the `.venv`-directory inside your build directory. It will also install required Python packages into this environment. You can see the current environment definition in the file `pyproject.toml` inside your build-directory.

To manually add Python packages run the following inside your build-directory:

```bash
poetry add python-package-name # e.g. poetry add numpy
```

To activate the environment run `poetry shell` (this opens a new shell) or for one-off commands run `poetry run command`, e.g. `poetry run pip list`.

### Poetry & Conan

You can also install Conan with Poetry (so you don't need to install it system-wide) with the CMake-option `OGS_USE_CONAN=auto`:

```bash
cmake ../ogs -DOGS_USE_CONAN=auto
```

### Poetry & Benchmarks

You can use the argument `PYTHON_PACKAGES` on `AddTest()` to specify additional Python package dependencies.

The following example would install the latest version of `numpy` and `pandas` version 0.1.2:

```cmake
AddTest(
    ...
    PYTHON_PACKAGES numpy pandas=0.1.2
    ...
)
```
