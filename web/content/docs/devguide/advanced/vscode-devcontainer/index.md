+++
date = "2018-02-26T11:00:13+01:00"
title = "Develop with VS Code Remote – Containers"
author = "Lars Bilke"
weight = 1061

[menu]
  [menu.devguide]
    parent = "advanced"
+++

[Visual Studio Code](https://code.visualstudio.com/) is a powerful text editor which can be expanded to a full featured integrated development environment (IDE) with plugins. With the  [Visual Studio Code Remote - Containers](https://code.visualstudio.com/docs/remote/containers) extension you can use prebuilt Docker container as the runtime environment for your code development. This page serves as an alternative to the [Getting Started]({{< relref "docs/devguide/getting-started/introduction.md" >}})-section. The following features make this a nice development environment:

- All ogs prerequisites (except of MPI / PETSc support, may come later) already setup.
- Code editing with auto-completion.
- Easy build configuration by simply selecting CMake presets.
- Start a debugging session with a click of a button.
- Hugo web server automatically running at [http://localhost:1313](http://localhost:1313) watching for changes in `web/content`.
- `zsh`-terminal with sane defaults (`oh-my-zsh`).
-
- Works also on a remote machine!

## Set-up prerequisites

- Install Git [as described in the Getting-Started section]({{< relref "prerequisites#step-install-git" >}}).
- Clone the OGS source code [as described in the Getting-Started section]({{< relref "get-the-source-code#option-clone-the-source-code-repository-with-git" >}}).
- Create a `build`-directory at the same level as the `ogs` source code directory.
- Install [Docker](https://docs.docker.com/get-docker/) (**or** have a Docker machine available where you have SSH access).
- Install [VS Code](https://code.visualstudio.com).
- Inside VS Code install the [Remote Development](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack) extension pack (includes *Remote – Containers* and *Remote – SSH* extensions).

## Open OGS source code inside the development container

- Open VS Code.
- Open the OGS source code folder.
- Press `F1` and type `reopen in [container]`, press `ENTER`.

This takes now some time as the container is downloaded from the registry. Once finished you should see the following in the bottom status bar:

![VS Code status bar](devcontainer_footer.png "VS Code status bar shows that we are connected the container `ogs-gcc-dev`.")

## Debug ogs

As an example use case we configure, build and debug the ogs executable.

CMake configuration is handled by using CMake presets which can be selected from the bottom status bar in VS Code. See the [Configure and build with CMake Presets in Visual Studio Code](https://github.com/microsoft/vscode-cmake-tools/blob/main/docs/cmake-presets.md) for details.

In the editor:

- Open `ogs.cpp`
- Set a [breakpoint](https://code.visualstudio.com/docs/editor/debugging#_breakpoints) in the first line of the `main()`-function (around line 58) by clicking on the left gutter in the editor window (a red dot marks the enabled breakpoint).

In the status bar:

1. Select the `debug` preset.
2. Click the `Build` button.
3. Select `ogs` as the launch target.
4. Click the bug icon.

After some seconds the status bar color changes orange to indicate an active debugging session:

![VS Code debugging session](devcontainer_debug.png "Active debugging session: Execution stopped at line 58 in `ogs.cpp`. Step over lines with the Debug toolbar buttons on top. See local variables on the left.")

## Use a remote machine

All of this works also via ssh. If you do not have Docker locally running but some server with Docker where you have ssh access to you can do the following:

- Open VS Code.
- Click the green button on the bottom left.
- Press `ENTER` (*Connect to host*).
- Choose the server.
- [Now you are in a [Remote – SSH](https://code.visualstudio.com/docs/remote/ssh)-session.]
- Open the OGS source code folder (you have to clone it before on that server).
- Press `F1` and type `reopen in [container]`, press `ENTER`.
