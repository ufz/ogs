+++
date = "2018-02-26T11:00:13+01:00"
title = "Executable Testing / Test Data"
author = "Lars Bilke"
weight = 1023

[menu]
  [menu.devguide]
    parent = "testing"
+++

## Introduction

An executable test can be run with several *wrappers*, e.g. with `valgrind`, `memcheck` or a simple `time`-measurement. Each wrapper run can then be validated (called *testers*), e.g. with file-`diff`s or -`grep`s. This can be configured in CMake:

```cmake
AddTest(
    NAME GroundWaterFlowProcess
    PATH Elliptic/quad_20x10_GroundWaterFlow
    EXECUTABLE ogs
    EXECUTABLE_ARGS quad_20x10_GroundWaterFlow.prj
    RUNTIME 35                                                        # optional
    WRAPPER time                                                      # optional
    TESTER diff                                                       # optional
    DIFF_DATA quad_20x10_constMat0.mesh.vtu quad_20x10_left_right.gml # optional
)
```

Tests are then run with `ninja ctest` or for more verbose output with `ctest -VV` (you may also use other [`ctest` options](https://cmake.org/cmake/help/v3.4/manual/ctest.1.html)). If the checker has some errors they are displayed. `RUNTIME` specifies the typical runtime in seconds on an Intel Xeon E5-2680 v2 @ 2.80 GHz with 500 GiB RAM (`envinf1`). Tests with a `RUNTIME > {{< dataFile "versions.ctest.large_runtime" >}}` are considered `LARGE`-tests.

The functionality is very flexible and more wrappers and checker can be added later on. e.g. for running some statistics on output files and comparing them with statistics from reference files.

<div class="note">

<h3>Run tests with CMake presets</h3>

Similar to the configure and build presets there are test presets, e.g. in your source-directory:

```bash
ctest --preset release                          # equivalent to running `ninja ctest` above
ctest --preset release -j 6 --label-regex Utils # run 6 tests in parallel which have a Utils label
```

**To sum up:** from a clean source directory you can fully configure, build and test OGS with these 3 commands:

```bash
cmake --preset release
cmake --build --preset release
ctest --preset release
```

</div>

## Test Data

Test data is stored in `Tests/Data`. Generated test output files should be found in `[build-dir]/Tests/Data`.

In the OGS-cli outputting to `[build-dir]/Tests/Data` is already handled (via the `-o` parameter). For other executables you have to implement this, e.g. a with parameter specifying the output directory.

In code `BaseLib::BuildInfo::data_path` (from `BuildInfo.h`) references the data source directory and `BaseLib::BuildInfo::data_binary_path` references the data output directory.

For adding new data files simply commit the new files as usual.

## Notebook testing

Full Jupyter Notebooks based workflows can be tested too. Create the notebook in `Tests/Data`. Configure output directory and try to use it for all outputs:

```python
import os

# On CI out_dir is set to the notebooks directory inside the build directory
# similar to regular benchmark tests. On local testing it will output to the
# notebooks source directory under a _out-subdirectory.
out_dir = os.environ.get('OGS_TESTRUNNER_OUT_DIR', '_out')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# ...
# Run ogs; get input data from current directory; write to `out_dir`
! ogs my_project.prj -o {out_dir} > {out_dir}/log.txt

# Verify results; on failure assert with:
assert False
# or
raise SystemExit()
```

Add new Python dependencies to `Test/Data/requirements.txt`.

### Register with CTest

Add the benchmark to CTest with e.g.:

```cmake
if(NOT OGS_USE_PETSC)
    NotebookTest(NOTEBOOKFILE Mechanics/Linear/SimpleMechanics.ipynb RUNTIME 10)
endif()
```

By registering notebooks [are automatically added]({{< relref "jupyter-docs" >}}) to the benchmark documentation page.

For local testing please note that you need to configure OGS with `OGS_USE_PIP=ON` (to automatically create a virtual environment in the build directory which is used by the notebook tests).

Then e.g. run all notebook test (`-R nb`) in parallel with:

```bash
source .venv/bin/activate # May need to be activated
ctest -R nb -j 4 --output-on-failure
```

<div class='note'>

### PyVista notebooks on headless Linux systems

PyVista (or VTK) requires a windowing environment for rendering. You can provide a virtual window with `xvfb-run`:

```bash
sudo apt install libgl1-mesa-glx xvbf # install xvfb

xvfb-run -a ctest [...] # provide a virtual window to the ctest-run
```

</div>
