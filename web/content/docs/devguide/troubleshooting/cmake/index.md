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

Also consider using the command line for CMake configuration as lots of CMake options (which modify requirements on third-party libraries) have to be set via the command line **before** CMake ran for the first time. E.g. when building with PETSc the following fails:

- Creating build directory
- Starting CMake Gui
- Click Configure
- Enable option `OGS_USE_PETSC`

When clicking *Configure*, CMake runs without `OGS_USE_PETSC` enabled and picks a VTK library which is not compiled for MPI usage. The picked VTK library is cached and never modified. So when enabling PETSc it still uses the old (and now wrong) VTK library.

So you have to use the CMake command line:

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

Check also [Conans troubleshooting page]({{< ref "conan.md" >}}) if you use Conan for dependencies.
