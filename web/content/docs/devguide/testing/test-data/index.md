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

## Performance tests

In order to prove the effectiveness of optimization features or to exclude
performance regressions it's good practice to implement performance tests. This can be
done with `pytest` and the [log parser](https://ogstools.opengeosys.org/stable/user-guide/logparser.html) from OGSTools.
In order to run such tests for OGS you need a set up environment.
See [this page](/docs/devguide/advanced/python-wheel) to learn how to do that.

One example performance test is the one [testing an optimization in the heat
transport BHE process assembly](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Python/test_bhe_assembly_matrix_cache_soil_elements.py?ref_type=heads).
It uses `pytest`'s abilities to modify environment variables and to capture the `stdout`
and `stderr` streams. It might serve as a starting point for your own performance
tests.

Performance tests relying on time measurements should be marked with the [`pytest` mark](https://docs.pytest.org/en/stable/how-to/mark.html) `performance_test` (it's a custom mark in the OGS test suite):

```py
@pytest.mark.performance_test()
def test_assembly_optimization():
    assert 1 == 0  # some meaningful test
```

That will ensure that this test is handled correctly in the OGS CI pipeline:
a failing performance tests will usually not be a strict failure in the OGS CI
pipelines, because on usual machines time measurements might be not meaningful
(e.g., if multiple concurrent processes run).
Therefore, in the OGS CI pipelines performance test failures will make the test
suite fail only on isolated CI agents, where no other jobs run at the same time.

That means: only mark tests relying on time measurements with
`@pytest.mark.performance_test()`!

### Performance test behaviour/selection

If a performance test failure is a strict failure can be controlled with the
environment variable `OGS_PERFORMANCE_TESTS_ALLOWED_TO_FAIL`. The behaviour is as follows:

* environment variable unset or set to `false`: performance test failure will be treated
  as a failure
* environment variable set to `true`: performance test failure will be treated
  as an expected failure (*xfail* in `pytest` terms) and will not make the test
  suite fail

Hence, the default is "treat as failure". On most CI agents the environment
variable is set to `true`, only on "isolated" CI machines it is set to `false`.

That means, that on you local machine by default failing performance tests will
make `pytest` fail.
Since your local machine is probably not an isolated machine you might want to
change that.
Thus, either you set the environment variable as described above or you deselect
performance tests entirely, e.g. (see also the [`pytest` documentation](https://docs.pytest.org/en/stable/example/markers.html#mark-examples)):

```sh
pytest -m 'not performance_test'
```

## Notebook tests

For notebook-based tests [see its dedicated page]({{< ref "jupyter-docs.md" >}}).
