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
- [pytest](https://docs.pytest.org/en/stable/ ) for testing Python functionality in OGS
- [nbconvert](https://nbconvert.readthedocs.io/en/latest/) for testing Jupyter Notebooks
- ...

Python packages are usually installed via `pip` inside an isolated environment (a virtual environment). We use a modern alternative to `pip`: [`uv`](https://docs.astral.sh/uv).

## Recommendation: Use the `OGS_USE_PIP=ON` CMake option when configuring OGS

Make sure that you installed the [uv](https://docs.astral.sh/uv/getting-started/installation/)-tool.

When configuring OGS with `OGS_USE_PIP=ON` `uv` creates a new virtual environment in the `.venv`-directory inside your build directory. It will also install required Python packages into this environment. For example, OGSTools or Jupyter will this way be available ready-made in a consistent configuration.

Example usage:

```bash
cmake --preset release -DOGS_USE_PIP=ON
```

When you want to use the python packages in this virtual environment, you need
to activate it by running

```bash
cd ../build/release          # the build directory
source .envrc                # Linux, Mac OR
Invoke-Expression .envrc.ps1 # Windows (PowerShell)
```

where `.envrc` is a script including `source .venv/bin/activate` for activating the virtual environment as well as setting some paths. If you have the [`direnv`](https://direnv.net)-tool installed and set up the virtual environment will be activated automatically upon changing into the build directory.

To manually add Python packages run the following inside your build-directory:

```bash
uv pip install python-package-name
```

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
