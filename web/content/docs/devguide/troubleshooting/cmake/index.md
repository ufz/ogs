+++
date = "2018-02-26T11:00:13+01:00"
title = "CMake"
author = "Lars Bilke"
weight = 1072

[menu]
  [menu.devguide]
    parent = "troubleshooting"
+++

## <i class="far fa-exclamation-triangle"></i> General advice

If something goes wrong when running CMake please try again with an **empty** or newly created build-directory! This is the very first thing you should try!

Please read the CMake output carefully. Often it will tell you what went wrong.

Also consider using the command line for CMake configuration as lots of CMake options (which modify requirements on third-party libraries) have to be set via the command line **before** CMake ran for the first time.

For example when building with PETSc the following fails:

- Creating build directory
- Starting CMake GUI
- Click Configure
- Enable option `OGS_USE_PETSC`

When clicking *Configure*, CMake runs without `OGS_USE_PETSC` enabled and picks a VTK library which is not compiled for MPI usage. The picked VTK library is cached and never modified. When enabling PETSc it still uses the old (and now wrong) VTK library.

Instead you should use the CMake command line:

```bash
mkdir build
cd build
cmake ../ogs -DOGS_USE_PETSC=ON
make
```

The following options are affected by this behavior and **should not be changed** after initially set (this means **do not** use `ccmake` or CMake-GUI for the first CMake run):

- `OGS_BUILD_GUI`
- `OGS_USE_PETSC`
- `CMAKE_BUILD_TYPE`
- `BUILD_SHARED_LIBS`

If the errors are related to TFEL, HDF5, PETSc or VTK, you may delete `~/.cache/CPM/_ext` and try again.

If there are incompatibilities to a system-installed library, e.g. HDF5, you can also use the additional cmake options, e.g. `-DOGS_BUILD_HDF5=ON`, to force a local library (in this case HDF5) build. These options are available:

- `OGS_BUILD_HDF5`
- `OGS_BUILD_TFEL`
- `OGS_BUILD_PETSC`
- `OGS_BUILD_VTK`

External dependencies are cached in `~/.cache/CPM/`.
If something goes wrong during the first build, then ogs will pick up the failed
build during subsequent cmake-runs without showing any error messages from the
dependency.
In that case you need to delete the corresponding folder in the CPM cache,
usually something like `~/.cache/CPM/_ext/TFEL/[some long hash]`.
Then the next cmake-run will retrigger a build of the dependency. This applies
to TFEL, HDF5, PETSc and VTK.
