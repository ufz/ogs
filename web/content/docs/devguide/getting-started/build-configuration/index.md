+++
date = "2018-02-23T15:28:13+01:00"
title = "Build configuration"
author = "Lars Bilke"
weight = 4

[menu]
  [menu.devguide]
    parent = "getting-started"
+++

## Overview

Before compiling the developer has to choose a configuration of the software. OGS comes in lots of different flavours (e.g. serial / parallelized), can be build with optional features or modules (e.g. MFront material models).

To separate source code from generated files such as compiled libraries, executables, test outputs and IDE projects we create build-directories. They can be placed arbitrarily. You can have as many build-directories as you like for e.g. different configurations but they will all use one source code directory. A typically directory structure:

- `ogs-source-code` (or simply `ogs`)
- `build` (should be placed outside the source directory)
  - `release`
  - `debug`

## Configure with CMake

For configuring a build the open source tool [CMake](http://www.cmake.org) is used. `CMakeLists.txt` files replace traditional Makefiles or IDE project files.

We provide CMake configuration presets defined in [CMakePresets.json](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/CMakePresets.json) for simple build configuration! Otherwise see [Configure manually](#option-configure-manually)). See the following table for commonly used presets:

### Available CMake presets

<!-- vale off -->

|               | Ninja[^1]     | Visual Studio    |
| ------------- | ------------- | ---------------- |
| CLI Release   | release       | msvc-release     |
| CLI Debug     | debug         | msvc-debug       |
| GUI Release   | release-gui   | msvc-release-gui |
| GUI Debug     | debug-gui     | msvc-debug-gui   |
| PETSc Release[^2] | release-petsc | -                |
| PETSc Debug[^2]   | debug-petsc   | -                |

<!-- vale on -->

[^1]: Requires the `ninja`-tool. See [install instructions]({{< ref "prerequisites.md#optional-install-ninja" >}}).
[^2]: Requires the `pkg-config`-tool. Can be installed via e.g. `apt` or `brew`.

### Configure with a preset

In the source directory run `cmake` with a preset:

```bash
# Usage: cmake --preset [preset-name]
cmake --preset release
```

This will create a build-directory outside the source tree (`../build/release`) with the default CMake options and the Release configuration.

Additionally you can pass any CMake variable or option with `-DVARIABLE_NAME=VALUE` (note the `-D` in front!) to the CMake command. You can also overwrite the generator with the `-G` parameter or the build-directory with the `-B` parameter (to see all available options just run `cmake --help`).

Also all the compiled files will be generated in this directory. This keeps the actual source code clean from intermediate files which are generated from the source code. Nothing inside the build directory will ever be version controlled because its contents can be regenerated anytime from the source code.

When you want to start over with a new configuration simply delete the build-directory, create a new one and reconfigure.

[See this]({{< ref "configuration-options" >}}) for a list of commonly used available options.

<div class='note'>

### User-defined presets

You can also create a `CMakeUserPresets.json` file in the root source directory with your own presets (this file is ignored by git):

```json
{
  "version": 1,
  "configurePresets": [
    {
      "name": "my-release",
      "inherits": "release",
      "cacheVariables": {
        "OGS_USE_PIP": "ON"
      }
    }
  ]
}

```

</div>

<div class='win'>

<div class='note'>

### Windows notes

#### <i class="far fa-check"></i> Ninja requirement: Use the Visual Studio command line

To use the Ninja build tool you need to configure in the Visual Studio command line. In the Start menu under *Visual Studio {{< dataFile "versions.minimum_version.msvc.year" >}}* you find a application link to *x64 Native Tools Command Prompt for VS {{< dataFile "versions.minimum_version.msvc.year" >}}*. This starts a command line setup for Visual Studio 64-bit. Now you can use a Ninja preset:

```bash
cmake --preset release
```

#### <i class="far fa-exclamation-triangle"></i> Multi-configuration with Visual Studio

OGS requires you to have one build directory which corresponds to one configuration. If you want to have e.g. a release and a debug build you need to create two build directories. This is nothing new to Linux / GCC user but differs to Visual Studios default behavior having just one build-folder / project with different configurations. A typical Visual Studio setup with both Release and Debug configurations would be initialized as follows:

```powershell
cmake --preset msvc-release # creates ../build/msvc-release/OGS.sln
cmake --preset msvc-debug   # creates ../build/msvc-debug/OGS.sln
```

Please also note that in Visual Studio you have to choose the correct configuration (i.e. when opening the solution-file in the release-folder you have to switch the Visual Studio configuration to **Release**)!

#### <i class="far fa-check"></i> Pro Tip: Use a better terminal application

Use the [Windows Terminal](https://apps.microsoft.com/detail/9N0DX20HK701) for a better terminal experience. It offers modern terminal features such as multiple tabs and panes.

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

<div class='note'>

### Using a python virtual environment

Additionally to

```bash
sudo apt-get install python3 python3-pip
```

(from [Set Up Prerequisites]({{< ref "prerequisites.md" >}})), do

```bash
sudo apt-get install python3-venv
```

and then you can use the `-DOGS_USE_PIP=ON` option:

```bash
cmake ../ogs -DOGS_USE_PIP=ON -DCMAKE_BUILD_TYPE="Release" -G Ninja
```

</div>

## Option: Configure with a visual tool

<div class='win'>

CMake comes with a graphical tool called `cmake-gui`. You can find it in the **Windows Start Menu**. First you need to set the source and build directory. Then click **Configure**. Now choose the generator to be used (e.g. **Visual Studio {{< dataFile "versions.minimum_version.msvc.number" >}} {{< dataFile "versions.minimum_version.msvc.year" >}}** for Visual Studio {{< dataFile "versions.minimum_version.msvc.year" >}}). Now choose your desired configuration options by toggling the corresponding checkboxes. Click **Configure** again. Click **Configure** often enough until the **Generate**-button becomes visible. Pressing **Generate** will finally generate the project files inside the chosen build directory.

</div>

<div class='linux'>

A more convenient way of running CMake on the command line is to use the `ccmake` tool. This is a shell tool but with some graphical user interface. If you want to use a preset configure with regular `cmake` first (`ccmake` **does not support presets**) and then run `ccmake` supplying the path to the build directory:

```bash
cmake --preset release
ccmake ../build/release
```

You are now presented the available configuration options. You can navigate in the list with the cursor keys and toggle / alter options with <kbd>Enter</kbd>. You may also press <kbd>T</kbd> to toggle (previously hidden) advanced options. Press <kbd>C</kbd> again until the **Generate**-option becomes visible. Press <kbd>G</kbd> to generate the project files and exit `ccmake`.

There is also the tool `cmake-gui` available (which supports presets), please see the Win-Tab for a description.

</div>

<div class='mac'>

Please see the Linux instructions!

</div>
