+++
date = "2018-02-26T11:00:13+01:00"
title = "Configuration options"
author = "Lars Bilke"
weight = 1033

[menu]
  [menu.devguide]
    parent = "advanced"
+++

## CMake options

Some of these options are enabled by default ("*Defaults* to *ON*") otherwise they must be explicitly set to *ON*.

### General

CMake switches to enable / disable parts of OGS.

- `OGS_BUILD_CLI` - Builds the simulator. *Defaults* to *ON*. If set to *OFF* all processes are also disabled.
- `OGS_BUILD_GUI` - Builds the Data Explorer. *Defaults* to *OFF*.
- `OGS_BUILD_TESTING` - Builds the test executables. *Defaults* to *ON*.
- `OGS_BUILD_UTILS` - Builds several utilities.
- `OGS_NO_EXTERNAL_LIBS` - Disables all external optional dependencies.
- `OGS_BUILD_PROCESS_X` - For enabling/disabling compilation of process `X`.
  Run the CMake-Gui to see a list of processes.
- `OGS_BUILD_PROCESSES` - A `;`-separated list specifying processes to build. *Defaults* to an *empty string*. This will alter the `OGS_BUILD_PROCESS_X`-options. For e.g. building just the two processes `HT` and `LIE`: `-DOGS_BUILD_PROCESSES="HT;LIE"`. Setting this variable back to an empty string **does not reset** the `OGS_BUILD_PROCESS_X`-options. You can also set it to *OFF* to disable all processes.

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
- `OGS_CPU_ARCHITECTURE` - Optimizes for the given CPU architecture see [-march](https://gcc.gnu.org/onlinedocs/gcc-4.5.3/gcc/i386-and-x86_002d64-Options.html)-flag. Defaults to `native`. For redistributable binaries set to `generic` on Linux and `core2` on Mac OS.
- `CMAKE_LIBRARY_SEARCH_PATH` - Additional library installation path, e.g. `/opt/local` or `C:/libs`
