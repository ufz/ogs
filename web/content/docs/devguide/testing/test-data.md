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
    RUNTIME 35                                                        # optional
    EXECUTABLE ogs                                                    # optional
    EXECUTABLE_ARGS quad_20x10_GroundWaterFlow.prj
    WRAPPER time                                                      # optional
    TESTER diff                                                       # optional
    DIFF_DATA quad_20x10_constMat0.mesh.vtu quad_20x10_left_right.gml # optional
)
```

Tests are then run with `make ctest` or for more verbose output with `ctest -VV` (you may also use other [ctest options](https://cmake.org/cmake/help/v3.4/manual/ctest.1.html)). If the checker has some errors they are displayed. `RUNTIME` specifies the typical runtime in seconds on `envinf1`. Tests with a `RUNTIME > 60` are considered `LARGE`-tests.

The functionality is very flexible and more wrappers and checker can be added later on. e.g. for running some statistics on output files and comparing them with statistics from reference files.

## Test Data

Test data is stored in `Tests/Data`. Generated test output files should be found in `[build-dir]/Tests/Data`.

In the OGS-cli outputting to `[build-dir]/Tests/Data` is already handled (via the `-o` parameter). For other executables you have to implement this, e.g. a with parameter specifying the output directory.

In code `BaseLib::BuildInfo::data_path` (from `BuildInfo.h`) references the data source directory and `BaseLib::BuildInfo::data_binary_path` references the data output directory.

For adding new data files simply commit the new files as usual.
