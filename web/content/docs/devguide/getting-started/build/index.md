+++
date = "2018-02-23T15:28:13+01:00"
title = "Build"
author = "Lars Bilke"
weight = 1005

[menu]
  [menu.devguide]
    parent = "getting-started"
+++

## Build the project

<div class='win'>

### Step: Build with Visual Studio

Open the OGS.sln (in the build directory) either by double-clicking it in the file browser or opening in Visual Studio via **File / Open / Project**.

On the project explorer right-click on **ogs** or **ogs-gui** and choose **Set as startup project**. Then press <kbd>F5</kbd> or click the green arrow to build and start the project.

#### About Visual Studio startup projects

The reason for this is that you can have only one sub-project of your Visual Studio Solution activated for debugging (e.g. by pressing <kbd>F5</kbd>). Per default this is the first sub-project (in the case of a CMake-generated project this is `ALL_BUILD`). This must be manually set to a sub-project which generates an executable (a library sub-project cannot be started). And because this setting is stored in user specific project file it cannot be set via CMake.

### How to work with CMake and Visual Studio

You can work normally in Visual Studio but remember that you have to make project changes in the `CMakeLists.txt`-file and not inside Visual Studio. So you can add a new source file within Visual Studios File menu but you have to add that file also to the CMake file. Every time you change a `CMakeLists.txt` and you build the project a new CMake run is automatically invoked. After that Visual Studio informs you that the project files were changed and it reloads them.
</div>

<div class='linux'>

### With ninja

To build with the `ninja` tool on the shell go to your previously configured build directory and type `ninja`:

```bash
cd ../build/release
ninja
```

<div class="note">

<h3><i class="fal fa-exclamation-triangle"></i> Exceeding memory consumption during compilation</h3>

By default, `ninja` uses maximum parallelism which may lead to a memory consumption exceeding the available memory (especially when compiling the processes).

Therefore, we recommend restricting the number of parallel compilations to be used to **one per 8 GB RAM**, e.g., if your machine has 16 GB RAM type `ninja -j 2` for a maximum of two parallel compilations.

Another option is to build just the process you are interested in by configuring with the CMake option `OGS_BUILD_PROCESSES`, e.g.: `cmake -DOGS_BUILD_PROCESSES=HT ....`. In this case you should be able to compile with full parallelization.

<h3>Build with CMake presets</h3>

Similar to the configure presets there are build presets, e.g. in your source-directory:

```bash
cmake --build --preset release             # equivalent to running ninja above
cmake --build --preset release -j 2 -t ogs # build the ogs-target with 2 threads
```

</div>

</div>

<div class='mac'>

See Linux-tab!

</div>

## Waiting

So now the build process is running... This can take some time because maybe there are external libraries which get automatically downloaded and compiled. This step is only done once per build directory, so subsequent builds will be much faster. See [this]({{< ref "cpm.md" >}}) for more.

## Finished

Congratulations you have finished the **Getting Started**-section!

Have a look at the other sections of this guide. Check out [Development Workflows]({{< ref "/docs/devguide/development-workflows/introduction.md" >}}) if you are interested in actively contributing to the project. The [Configuration Options]({{< ref "configuration-options" >}})-page shows you all available build customizations. Go ahead!
