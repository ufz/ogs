+++
date = "2021-07-19T11:34"
title = "MFront"
author = "Lars Bilke"
weight = 1043

[menu]
  [menu.devguide]
    parent = "packages"
+++

MFront / TFEL support is enabled by the CMake-option `-DOGS_USE_MFRONT=ON` and can be installed

- system-wide,
- in the user directory or
- automatically inside the ogs build-directory.

## System-wide install

If you install it system-wide it will get picked up by CMake (when the `mfront`-executable is in the `PATH`). Make sure to install a compatible version! Currently OGS requires the TFEL branch [rliv-{{< dataFile "versions.minimum_version.tfel-rliv" >}}](https://github.com/thelfer/tfel/tree/rliv-{{< dataFile "versions.minimum_version.tfel-rliv" >}}). Check the CMake output for version information of MFront / TFEL.

## Installation in the user directory

If you need another version of MFront / TFEL (e.g. for model development) as required for OGS one solution is to install MFront / TFEL in your user-directory with a version-suffix, see the [official documentation](http://tfel.sourceforge.net/install.html#sec:QuickUbuntu) for instructions! In this case you can use the mfront executable by using the binary name with the suffix, e.g. `mfront-4.0.0-dev` and OGS will automatically download and build a compatible MFront / TFEL version inside the build directory.

## Automatic installation

If you do not want to care about all this simply do **not** install MFront / TFEL and OGS will automatically download and build a compatible MFront / TFEL version inside the build directory.
