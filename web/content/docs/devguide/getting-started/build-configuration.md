+++
date = "2018-02-23T15:28:13+01:00"
title = "Build configuration"
author = "Lars Bilke"
weight = 1004

[menu]
  [menu.devguide]
    parent = "getting-started"
+++

## Overview

Before compiling the developer has to choose a configuration of the software. OGS comes in lots of different flavours (e.g. serial / parallelized), can be build with optional features (e.g. with Python scripting support) or modules (e.g. MFront material models).

To separate source code from generated files such as compiled libraries, executables, test outputs and IDE projects we create build-directories. They can be placed arbitrarily. You can have as many build-directories as you like for e.g. different configurations but they will all use one source code directory. A typically directory structure:

- `ogs-source-code` (or simply `ogs`)
- `build` (should be placed outside the source directory)
  - `release`
  - `debug`

## Configure with CMake

For configuring a build the open source tool [CMake](http://www.cmake.org) is used. `CMakeLists.txt` files replace traditional Makefiles or IDE project files.

We provide CMake configuration presets defined in [CMakePresets.json](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/CMakePresets.json) for simple build configuration (**Note:** Requires CMake $\geq$ {{< dataFile "versions.tested_version.cmake" >}}! Otherwise see [Configure manually](#option-configure-manually)). See the following table for commonly used presets:

### Available CMake presets

|               | Ninja[^1]     | Visual Studio    |
| ------------- | ------------- | ---------------- |
| CLI Release   | release       | msvc-release     |
| CLI Debug     | debug         | msvc-debug       |
| GUI Release   | release-gui   | msvc-release-gui |
| GUI Debug     | debug-gui     | msvc-debug-gui   |
| PETSc Release | release-petsc | -                |
| PETSc Debug   | debug-petsc   | -                |

[^1]: Requires the `ninja`-tool. See [install instructions]({{< ref "prerequisites.md#optional-install-ninja" >}}).

### Configure with a preset

In the source directory run `cmake`:

```bash
# Usage: cmake -S [path-to-source] --preset=[preset-name]
cmake -S . --preset release
```

This will create a build-directory outside the source tree (`../build/release`) with the default CMake options and the Release configuration.

Additionally you can pass any CMake variable or option with `-DVARIABLE_NAME=VALUE` (note the `-D` in front!) to the CMake command. You can also overwrite the generator with the `-G` parameter or the build-directory with the `-B` parameter (to see all available options just run `cmake --help`)

Also all the compiled files will be generated in this directory. This keeps the actual source code clean from intermediate files which are generated from the source code. Nothing inside the build directory will ever be version controlled because its contents can be regenerated anytime from the source code.

When you want to start over with a new configuration simply delete the build-directory, create a new one and reconfigure.

[See this]({{< ref "configuration-options" >}}) for a list of commonly used available options.

<div class='note'>

### User-defined presets

You can create a `CMakeUserPresets.json` file in the root source directory with your own presets (this file is ignored by git):

```json
{
  "version": 1,
  "configurePresets": [
    {
      "name": "my-release",
      "inherits": "release",
      "cacheVariables": {
        "OGS_USE_POETRY": "OFF"
      }
    }
  ]
}

```

</div>

<div class='win'>

<div class='note'>

### Windows notes:

#### <i class="far fa-check"></i> Ninja requirement: Use the Visual Studio command line

To use the Ninja build tool you need to configure in the Visual Studio command line. In the Start menu under *Visual Studio {{< dataFile "versions.minimum_version.msvc.year" >}}* you find a application link to *x64 Native Tools Command Prompt for VS {{< dataFile "versions.minimum_version.msvc.year" >}}*. This starts a command line setup for Visual Studio 64-bit. Now you can use a Ninja preset:

```bash
cmake -S . --preset=release
```

#### <i class="far fa-exclamation-triangle"></i> Multi-configuration with Conan and Visual Studio

With Conan one build directory corresponds to one configuration. If you want to have e.g. a release and a debug build you need to create two build directories. This is nothing new to Linux / GCC user but differs to Visual Studios default behavior having just one build-folder / project with different configurations. A typical Visual Studio setup with both Release and Debug configs would be initialized as follows:

```bash
cmake -S . --preset=msvc-release
cmake -S . --preset=msvc-debug
```

Please also note that in Visual Studio you have to choose the correct configuration (i.e. when opening the solution-file in the release-folder you have to switch the Visual Studio configuration to **Release**)!

#### <i class="far fa-check"></i> Pro Tip: Use a better terminal application

Use the [Windows Terminal](https://www.microsoft.com/en-us/p/windows-terminal/9n0dx20hk701) for a better terminal experience. It offers modern terminal features such as multiple tabs and panes.

</div>

</div>

## Option: Configure manually

If you cannot use CMake presets (e.g. when your CMake installation does not support it) manually create a build directory and run CMake from within with all required parameters, e.g:

```bash
# in ogs source-directory
mkdir -p ../build/release
cd ../build/release
cmake ../../ogs -G Ninja -DCMAKE_BUILD_TYPE=Release
```

<div class='note'>

### Using a different compiler

Set the `CC` and `CXX` environment variables, e.g.:

```bash
CC=mpicc CXX=mpic++ cmake ../ogs -G Ninja -DCMAKE_BUILD_TYPE=Release -DOGS_USE_PETSC=ON
```

</div>

## Option: Configure with a visual tool

<div class='win'>

CMake comes with a graphical tool called **cmake-gui**. You can find it in the **Windows Start Menu**. First you need to set the source and build directory. Then click **Configure**. Now choose the generator to be used (e.g. **Visual Studio {{< dataFile "versions.minimum_version.msvc.number" >}} {{< dataFile "versions.minimum_version.msvc.year" >}}** for Visual Studio {{< dataFile "versions.minimum_version.msvc.year" >}}). Now choose your desired configuration options by toggling the corresponding checkboxes. Click **Configure** again. Click **Configure** often enough until the **Generate**-button becomes visible. Pressing **Generate** will finally generate the project files inside the chosen build directory.

</div>

<div class='linux'>

A more convenient way of running cmake on the command line is to use the `ccmake` tool. This is a shell tool but with some graphical user interface. To use it just run `ccmake` instead of `cmake`:

```bash
ccmake -S . --preset=release
```

First press <kbd>C</kbd> to **Configure**. You are now presented the available configuration options. You can navigate in the list with the cursor keys and toggle / alter options with <kbd>Enter</kbd>. You may also press <kbd>T</kbd> to toggle (previously hidden) advanced options. Press <kbd>C</kbd> again until the **Generate**-option becomes visible. Press <kbd>G</kbd> to generate the project files and exit `ccmake`.

There is also the tool `cmake-gui` available, please see the Win-Tab for a description.

</div>

<div class='mac'>

Please see the Linux instructions!

</div>
