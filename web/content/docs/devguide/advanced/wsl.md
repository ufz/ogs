+++
date = "2020-05-15T10:46"
title = "Windows Subsystem for Linux"
author = "Lars Bilke"
weight = 1042

[menu]
  [menu.devguide]
    parent = "advanced"
+++

The [Windows Subsystem for Linux](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux) is an alternative way to setup a complete development environment on Windows. It offers a Linux environment with a seamless bridge to the Windows world. We recommend this setup for Windows developers.

## Setup

- Install WSL by following [this guide](https://www.windowscentral.com/install-windows-subsystem-linux-windows-10). **Important:** Choose **Ubuntu 20.04 LTS** as the Linux distribution. Other distros may not have a sufficient compiler.
- Follow the [developer guide for Linux]({{< ref "prerequisites.md" >}}) from now on.

## Using Visual Studio Code as IDE

You can use the native Windows Visual Studio Code IDE (VS Code) for developing in the WSL. It offers code completion, CMake integration, an integrated debugger, git integration and more.

### Setup

- In the WSL shell run `sudo apt update && sudo apt install -y ninja-build gdb`.
- On Windows [install Visual Studio Code](https://code.visualstudio.com/docs/setup/windows)
- Open VS Code and install the [VS Code Remote Development](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack) extension
- In the WSL shell go to your OGS source code and run `code .`. This will open the source code in VS Code.
- Install the [C/C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools), [CMake](https://marketplace.visualstudio.com/items?itemName=twxs.cmake) and [CMake Tools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools) extensions

### Build and debug OGS

- Configure your project by following [this guide](https://vector-of-bool.github.io/docs/vscode-cmake-tools/getting_started.html#configuring-your-project). Select `GCC 9.3.0` as the [CMake kit](https://vector-of-bool.github.io/docs/vscode-cmake-tools/kits.html#kits).
- Follow [this guide](https://vector-of-bool.github.io/docs/vscode-cmake-tools/debugging.html#selecting-a-launch-target) for debug a target. Select `ogs` as the debug target.

## Additional notes

The filesystem of the WSL is not inside your regular user directories. You can find it by running `explorer.exe .` inside the WSL shell. It should be something like `\\wsl$\Ubuntu-20.04\home\[username]\...`.

You can also run OGS inside the WSL with benchmarks located in your regular Windows directories. You regular filesystem can be accessed inside WSL with the `/mnt/c/`-prefix. E.g. to run an OGS benchmark:

```bash
bin/ogs -o _out /mnt/c/Users/[username]/ogs-src/Tests/Data/Elliptic/square_1x1_SteadyStateDiffusion/square_1e0.prj
```
