+++
date = "2022-02-09T11:00:13+01:00"
title = "Python wheel development"
author = "Lars Bilke"
weight = 1068

[menu]
  [menu.devguide]
    parent = "advanced"
+++

## Local setup

Python wheel builds are driven by [scikit-build](https://scikit-build.readthedocs.io/en/latest/) which basically is a `setuptools`-wrapper for CMake-based projects.

The entrypoint is `setup.py` in the root directory. It uses the `wheel` CMake preset (or `wheel-win` on Windows).

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
```

The `pip install`-step starts a new CMake-based ogs build in `_skbuild`-subdirectory (inside the source code) using the `wheel`-preset. When the CMake build is done it installs the wheel into the activated virtual environment and you can interact with it, e.g.:

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

If you make modifications you need to run `pip install .[test]` again (or for temporary modifications you can directly edit inside the virtual environment, e.g. in `.venv/lib/python3.10/site-packages/ogs`).

The contents of `_skbuild/[platform-specific]/cmake-install` will make up the wheel.

## CI

For generating the various wheels for different Python versions and platforms [`cibuildwheel`](https://cibuildwheel.readthedocs.io/en/stable/) is used.

You can test it locally with, e.g. only building for Python 3.10:

```bash
CIBW_BUILD="cp310*" pipx run cibuildwheel
```

Please note that on Linux `cibuildwheel` runs the builds inside [manylinux](https://github.com/pypa/manylinux) Docker containers. On other platforms the build happens with native tools. See the [cibuildwheel docs](https://cibuildwheel.readthedocs.io/en/stable/#how-it-works) for more information.

Wheels are generated in the `wheelhouse/`-folder.

`cibuildwheel` is configured in `pyproject.toml`:

```toml
[tool.cibuildwheel]
...
```
