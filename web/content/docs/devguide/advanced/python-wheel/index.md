+++
date = "2022-02-09T11:00:13+01:00"
title = "Python bindings development"
author = "Lars Bilke"
weight = 1068

[menu]
  [menu.devguide]
    parent = "advanced"
+++

There are two ways to build OGS with Python bindings:

- A regular build with `OGS_USE_PIP=ON`
- A wheel build

## Regular build

- Configure your build with `-DOGS_USE_PIP=ON`
- Build

To run the Python-based tests:

- Source `build-dir/.envrc` (i.e., run `source build-dir/.envrc` in your terminal. This activates the virtual environment and sets `OGS_USE_PATH` environment variable).
- Run e.g. `pytest` from inside the source directory. On Windows run `Invoke-Expression build-dir/.envrc.ps1` in your PowerShell.
- Run e.g. `pytest` from inside the source directory or `pytest ../path/to/source/dir`.
- Alternatively you can also run the Python-based with `ctest`, e.g. `ctest -R pytest`.

If you make modifications on the C++ side you need to run `make` or `ninja` in the build directory again. If you modify Python code you need to run `cmake .` again in the build directory. Modifications on the Python tests are immediately available to `pytest`.

To get the output of a specific test:

```bash
pytest --capture=tee-sys ./Tests/Python/test_simulator_mesh_interface.py
```

## Wheel build

Python wheel builds are driven by [scikit-build-core](https://scikit-build-core.readthedocs.io) which basically is a `setuptools`-wrapper for CMake-based projects.

The entry point is `pyproject.toml` in the root directory. It uses the `wheel` CMake preset. The preset can be overridden and even other CMake options can be passed via the environment variable `CMAKE_ARGS` (see also the [scikit-build-core documentation](https://scikit-build-core.readthedocs.io/en/latest/configuration/index.html#configuring-cmake-arguments-and-defines)).

You can locally develop and test with the following setup:

```bash
# Create a virtual environment inside your source directory
python3 -m venv .venv
# Activate the environment
source .venv/bin/activate
# Install (build) the local Python project
pip install -v .[test]
...
Successfully installed ogs-6.4.2.dev1207

# To build with additional CMake arguments, e.g.:
pip install -v .[test] --config-settings=cmake.define.OGS_BUILD_PROCESSES=SteadyStateDiffusion"
```

The `pip install`-step starts a new CMake-based ogs build in `_skbuild`-subdirectory (inside the source code) using the `wheel`-preset. When the CMake build is done it installs the wheel into the activated virtual environment and you can interact with it.

The contents of `_skbuild/[platform-specific]/cmake-install` will make up the wheel.

### Testing

```bash
# Run python tests
pytest
============================================== test session starts ===============================================
platform darwin -- Python 3.10.6, pytest-7.1.3, pluggy-1.0.0
rootdir: ~/code/ogs/ogs, configfile: pyproject.toml, testpaths: Tests
collected 2 items

Tests/Python/test_cli.py .                                                                                 [ 50%]
Tests/Python/test_simlator.py .                                                                            [100%]

=============================================== 2 passed in 0.55s ================================================

# Start the python interpreter
python3
>>> import ogs.simulator as sim
>>> sim.initialize(["", "--help"])
```

If you make modifications on the C++ side you need to run `pip install .[test]` again. Modifications on the Python tests are immediately available to `pytest`.

To get the output of a specific test:

```bash
pytest --capture=tee-sys ./Tests/Python/test_simulator_mesh_interface.py
```

### Module structure

Python modules added in CMake via `pybind11_add_module()` should only contain the binding code and other (helper code) which is used by that module only! If you need helper code which is used by several modules (e.g. `OGSMesh`-class used in `mesh`- and `simulator`-module) it needs to be defined in a regular library and linked to the modules.

If you don't do this you will get unresolved externals between the modules and when you try to link them together you will get runtime errors or it doesn't link at all (at least on Windows, because of the `.pyd` library format).

## CI

For generating the various wheels for different Python versions and platforms [`cibuildwheel`](https://cibuildwheel.readthedocs.io/en/stable/) is used.

You can test it locally with, e.g. only building for Python 3.10:

```bash
CIBW_BUILD="cp310*" pipx run cibuildwheel
```

Please note that on Linux `cibuildwheel` runs the builds inside [`manylinux`](https://github.com/pypa/manylinux) Docker containers. On other platforms the build happens with native tools. See the [`cibuildwheel` docs](https://cibuildwheel.readthedocs.io/en/stable/#how-it-works) for more information.

Wheels are generated in the `wheelhouse/`-folder.

`cibuildwheel` is configured in `pyproject.toml`:

```toml
[tool.cibuildwheel]
...
```
