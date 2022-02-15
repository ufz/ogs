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

### General

CMake switches to enable / disable parts of OGS.

- `OGS_BUILD_CLI` - Builds the simulator. *Defaults* to *ON*. If set to *OFF* all processes are also disabled (see variable `OGS_BUILD_PROCESSES` below).
- `OGS_BUILD_GUI` - Builds the Data Explorer. *Defaults* to *OFF*.
- `OGS_BUILD_TESTING` - Builds the test executables. *Defaults* to *ON*.
- `OGS_BUILD_UTILS` - Builds several utilities.
- `OGS_BUILD_PROCESS_X` - For enabling/disabling compilation of process `X`.
  Run the CMake-Gui / ccmake to see a list of processes.
- `OGS_BUILD_PROCESSES` - A `;`-separated list specifying processes to build, e.g. `-DOGS_BUILD_PROCESSES="HT;LIE"`. Can be set to *ON* which means all processes are built or can be set to *OFF* to disable all processes. **Attention:** Setting this variable overrides individual `OGS_BUILD_PROCESS_X`-variables! This option is mainly used for CI and automation. Also the value of this variable is not cached.

### Debugging

- `CMAKE_BUILD_TYPE` - Defaults to `Debug` which builds with debugging infos, set to `Release` for an optimized build.
- `OGS_PROFILE` - Builds with profiling flags (`-pg`).
- `OGS_CMAKE_DEBUG` - Prints out the values of all defined CMake variables at CMake configuration time.

### Optimization

- `CMAKE_BUILD_TYPE` - Set to `Release` to build with optimization flags, set to `Debug` for debugging.

### Testing

- `OGS_COVERAGE` - Enables code coverage measurements with gcov/lcov. TODO

### Advanced options

- `OGS_CXX_FLAGS` - Appends user-given compiler flags. Note that existing (CMake-given) flags are not replaced.
- `OGS_PACKAGE_ADDITIONAL_BINARIES` - Package additional binaries (given as a `;`-separated list with relative paths to `CMAKE_BINARY_DIR`) into redistributables. Is used for bundling the OGS File Converter with the Data Explorer.
- `OGS_CPU_ARCHITECTURE` - Optimizes for the given CPU architecture see [-march](https://gcc.gnu.org/onlinedocs/gcc-4.5.3/gcc/i386-and-x86_002d64-Options.html)-flag. Defaults to `native`. For redistributable binaries set to `generic` on Linux and `core2` on Mac OS. Can be disabled when set to `OFF`.
- `CMAKE_LIBRARY_SEARCH_PATH` - Additional library installation path, e.g. `/opt/local` or `C:/libs`
- `OGS_DEPENDENCY_VERSIONS` - Overwrite individual entries in `web/data/versions.json`. Should be quoted and `;`-separated, e.g.: `-DOGS_DEPENDENCY_VERSIONS="minimum_version.petsc=3.16.2;ctest.large_runtime=120"`.
