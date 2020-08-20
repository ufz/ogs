+++
date = "2018-02-26T11:00:13+01:00"
title = "Build with ninja"
author = "Lars Bilke"
weight = 1031

[menu]
  [menu.devguide]
    parent = "advanced"
+++

## Install ninja

<div class='win'>

Download *ninja.zip* from the [latest GitHub release](https://github.com/ninja-build/ninja/releases/latest). Unzip it and make sure that the directory containing the `ninja.exe` is the `PATH`-environment variable!

</div>

<div class='linux'>

Install with your package manager:

```bash
sudo apt-get install ninja-build
```

</div>

<div class='mac'>

Install via Homebrew:

```bash
brew install ninja
```

</div>

## Configure for ninja and build

Run CMake with the ninja-generator:

```bash
cmake ../source/dir -G Ninja
ninja
```

<div class='note'>

### <i class="far fa-exclamation-triangle"></i> Visual Studio remarks

When you configure with the Ninja generator you have to run CMake from the appropriate Visual Studio Command Prompt! From there you can both use `cmake` as well as `cmake-gui` which starts the GUI. In the GUI select the `Ninja` generator and leave the toggle `Use default native compilers` on.

</div>
