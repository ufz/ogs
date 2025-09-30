+++
date = "2018-02-26T11:00:13+01:00"
title = "Configuration options"
author = "Lars Bilke"
weight = 1060

[menu]
  [menu.devguide]
    parent = "advanced"
+++

## CMake options

Some of these options are enabled by default ("*Defaults* to *ON*") otherwise they must be explicitly set to *ON*.

You can pass any CMake variable or option with `-DVARIABLE_NAME=VALUE` (note the
`-D` in front) to the CMake command.
You can also overwrite the generator with the `-G` parameter or the
build-directory with the `-B` parameter.

### General

CMake switches to enable / disable parts of OGS.

- `OGS_BUILD_CLI` - Builds the simulator. *Defaults* to *ON*. If set to *OFF* all processes are also disabled (see variable `OGS_BUILD_PROCESSES` below).
- `OGS_BUILD_GUI` - Builds the Data Explorer. *Defaults* to *OFF*.
- `OGS_BUILD_TESTING` - Builds the test executables. *Defaults* to *ON*.
- `OGS_BUILD_UTILS` - Builds several utilities.
- `OGS_BUILD_PROCESS_X` - For enabling/disabling compilation of process `X`.
  Run the CMake-GUI / `ccmake` to see a list of processes.
- `OGS_BUILD_PROCESSES` - A `;`-separated list specifying processes to build, e.g. `-DOGS_BUILD_PROCESSES="HT;LIE"`. Can be set to *ON* which means all processes are built or can be set to *OFF* to disable all processes. **Attention:** Setting this variable overrides individual `OGS_BUILD_PROCESS_X`-variables! This option is mainly used for CI and automation. Also the value of this variable is not cached.

### Python virtual environment

- `OGS_USE_PIP` - Creates a python virtual environment in the .venv-directory
  inside the build directory, that includes useful packages like ogstools,
  Jupyter and more. Defaults to *OFF*. See also
  https://www.opengeosys.org/docs/devguide/packages/python-env/.

### Debugging

- `CMAKE_BUILD_TYPE` - Defaults to `Debug` which builds with debugging information, set to `Release` for an optimized build.
- `OGS_PROFILE` - Builds with profiling flags (`-pg`).
- `OGS_CMAKE_DEBUG` - Prints out the values of all defined CMake variables at CMake configuration time.

### Optimization

- `CMAKE_BUILD_TYPE` - Set to `Release` to build with optimization flags, set to `Debug` for debugging.

### Testing

- `OGS_COVERAGE` - Enables code coverage measurements with `gcov` / `lcov`. TODO

### Advanced options

- `OGS_CXX_FLAGS` - Appends user-given compiler flags. Note that existing (CMake-given) flags are not replaced.
- `OGS_PACKAGE_ADDITIONAL_BINARIES` - Package additional binaries (given as a `;`-separated list with relative paths to `CMAKE_BINARY_DIR`) into redistributables. Is used for bundling the OGS File Converter with the Data Explorer.
- `OGS_CPU_ARCHITECTURE` - Optimizes for the given CPU architecture see [-march](https://gcc.gnu.org/onlinedocs/gcc-4.5.3/gcc/i386-and-x86_002d64-Options.html)-flag. Defaults to `native`. For redistributable binaries set to `generic` on Linux and `core2` on Mac OS. Can be disabled when set to `OFF`.
- `CMAKE_LIBRARY_SEARCH_PATH` - Additional library installation path, e.g. `/opt/local` or `C:/libs`
- `OGS_DEPENDENCY_VERSIONS` - Overwrite individual entries in `web/data/versions.json`. Should be quoted and `;`-separated, e.g.: `-DOGS_DEPENDENCY_VERSIONS="minimum_version.petsc=3.16.2;ctest.large_runtime=120"`.
- `OGS_USE_MKL` - Enables MKL support. Requires MKL (as part of the Intel oneAPI toolkit) to be [installed](https://www.intel.com/content/www/us/en/developer/articles/guide/installation-guide-for-oneapi-toolkits.html) on the system.

  Before configuring with CMake you also need to source the `setvars.sh`-script. On a typical installation run this: `source /opt/intel/oneapi/setvars.sh` (Windows: `& "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"`). This will set the `MKLROOT` environment variable and also the `DYLD_LIBRARY_PATH` to contain the MKL library locations (both variables will be checked by CMake).

  To enable 64-bit array indices in MKL add `-DMKL_INTERFACE=ilp64` on the first CMake run (with a clean CMake cache) but this seems [not supported by Eigen](https://libeigen.gitlab.io/docs/TopicUsingIntelMKL.html).

- `OGS_EIGEN_PARALLEL_BACKEND` - Defaults to `OpenMP`. Defaults to `MKL` when `OGS_USE_MKL=ON`. May be set to `OpenMP` when MKL is on which also enables OpenMP multithreading. This pulls in another OpenMP implementation besides the Intel MKL OpenMP which is an **experimental feature!**. If you want to use OpenMP parallelized assembly with MKL enabled you also need to explicitly set `OGS_EIGEN_PARALLEL_BACKEND=OpenMP`!
- `OGS_PETSC_CONFIG_OPTIONS` – Can be used to build a specific PETSc configuration, arguments are passed to PETSc's `./configure`-command. E.g. to build with MUMPS:

  ```bash
  cmake --preset release-petsc -DOGS_PETSC_CONFIG_OPTIONS="--download-superlu --download-superlu_dist --download-scalapack --download-mumps --with-cc=mpicc --with-cxx=mpic++ --with-fc=mpif90"
  ```
