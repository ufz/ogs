+++
date = "2025-05-11"
title = "Debugging OpenGeoSys with VS Code"
author = "Mostafa Mollaali, Lars Bilke"
weight = 1062

[menu]
  [menu.devguide]
    parent = "advanced"
+++



[Visual Studio Code](https://code.visualstudio.com/) is a powerful text editor that can be extended into a full-featured IDE using plugins. This guide walks you through setting up debugging for OpenGeoSys (OGS) on **macOS** and **Linux**.

---

## 1. Prerequisites

### Install system tools

**On macOS:**

```bash
xcode-select --install
````

**On Linux:**

```bash
sudo apt update
sudo apt install build-essential cmake gdb lldb ninja-build
```

### Clone and build OGS in debug mode

```bash
git clone https://gitlab.opengeosys.org/ogs/ogs.git
cd ogs
cmake --preset debug
ninja -C build/debug ogs
```

This creates a debug build in `build/debug/`.

---

## 2. VS Code Setup

### Install VS Code

Download from [https://code.visualstudio.com](https://code.visualstudio.com)

### Install required extensions

In the Extensions view (`Ctrl+Shift+X` on Linux, `Cmd+Shift+X` on macOS), install:

* [C/C++ Extension Pack](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools-extension-pack)
* [CMake Tools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools)
* [CodeLLDB](https://marketplace.visualstudio.com/items?itemName=vadimcn.vscode-lldb)
  *(required on macOS, optional on Linux if using LLDB)*

---

## 3. Project configuration

### Open the project

Open the  source code `ogs/` folder in VS Code.

### Optional: Add a build task

To build OGS from within VS Code:

```bash
mkdir -p .vscode
cat > .vscode/tasks.json <<EOF
{
  "version": "2.0.0",
  "tasks": [
    {
      "label": "Build OGS (CMake debug)",
      "type": "shell",
      "command": "cmake --preset debug && ninja -C ../build/debug ogs",
      "group": {
        "kind": "build",
        "isDefault": true
      }
    }
  ]
}
EOF
```

To run the build task, press `Shift + Command + B` on macOS or `Ctrl+Shift+B` on Linux.

---

### Configure debugging

Create `.vscode/launch.json`:

```bash
[ "${PWD##*/}" = ".vscode" ] || cd .vscode
cat > launch.json <<EOF
{
   "version": "0.2.0",
  "configurations": [
    {
      "name": "Debug with LLDB",
      "type": "lldb",
      "request": "launch",
      "program": "${workspaceFolder}/../build/debug/bin/ogs",
      "args": [
        "${workspaceFolder}/Tests/Data/LIE/Mechanics/coulomb_load_path.prj"
      ],
      "cwd": "${workspaceFolder}/../build/debug",
      "stopOnEntry": false,
      "externalConsole": false
    }
  ]
}
EOF
```

---

## 4. Start a debug session

1. Open the file `Applications/CLI/ogs.cpp`.
2. Set a breakpoint at the first line inside the `main()` function.
![Start debugging](vscode-debugging-1.png)
3. Start debugging by clicking the green "Run and Debug" icon, or press `F5` (Linux) / `Fn + F5` (macOS).
4. Choose the debug configuration by its name set in `launch.json` (e.g., `Debug with LLDB`).
![Running debug session](vscode-debugging-3.png)
5. The debugger will stop at your breakpoint. You can now step through code, inspect variables, and use the debug toolbar to control execution.
![Running debug session](vscode-debugging-2.png)
