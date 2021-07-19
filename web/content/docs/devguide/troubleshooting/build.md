+++
date = "2018-02-26T11:00:13+01:00"
title = "Build"
author = "Lars Bilke"
weight = 1073

[menu]
  [menu.devguide]
    parent = "troubleshooting"
    identifier = "build-troubleshooting"
+++

## Visual Studio out-of-heap or stackoverflow errors

<div class='note'>

**Note:** To prevent this you can also use the [WSL setup]({{< ref "wsl.md" >}}).

</div>

The compilation especially of the processes in Release-config can be very memory hungry. Using dynamic Eigen shape matrices can reduce memory usage:

```bash
cmake . -DOGS_EIGEN_DYNAMIC_SHAPE_MATRICES=ON
```

You should also have at least 8 GB of RAM and even this can be not enough when compiling on multiple cores (which is the default). To build on only one core run the following in your build-directory:

```bash
cmake --build . --config Release -j 1
```

If this still fails you can disable building of the failing processes, e.g.:

```bash
cmake . -DOGS_BUILD_PROCESS_HT=OFF
cmake --build . --config Release -j 1
```

## Compiler crash on Linux / macOS

The same as described above for Windows can happen on other platforms if there is not enough RAM to compile in parallel. See the note on [Exceeding memory consumption during compilation]({{< ref "docs/devguide/getting-started/build.md#build-the-project" >}}).
