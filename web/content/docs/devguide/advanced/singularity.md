+++
date = "2018-09-21T11:00:13+01:00"
title = "Singularity"
author = "Lars Bilke"
weight = 1036

[menu]
  [menu.devguide]
    parent = "advanced"
+++

<div class='note'>

### Important note

This page describes how to **build** OGS with the help of a Linux container (for **developers**). To **run** OGS with the help of a container go to the [user guide]({{< ref "container.md" >}}) (for **developers**).

</div>

## Introduction

[Singularity](https://www.sylabs.io) is a Linux container runtime similar to Docker. Key advantages over Docker are

- Container don't run with root privileges
- You are the same user with the same privileges inside the container as on the host
- Container can run on HPC systems and seamlessly integrate with resource managers and MPI
- Container can leverage NVidia GPUs

Singularity per default mounts your home directory and also passes your current working directory when starting a container. Therefore it is easy to use it for development.

### Prerequisites

- Linux and [Mac](https://sylabs.io/singularity-desktop-macos/) only!
  - Note: You can run Singularity **on Windows** inside a lightweight virtual machine using [WSL 2](https://docs.microsoft.com/en-us/windows/wsl/install-win10). Under Ubuntu use Homebrew package `singularity`, under CentOS install via [EPEL](https://sylabs.io/guides/3.0/user-guide/installation.html#install-the-centos-rhel-package-using-yum).
- [Install Git]({{< ref "prerequisites" >}}/#step-install-git)
- [Install Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps)

### Build OGS inside a container

```bash
[git clone ogs]
singularity pull docker://registry.opengeosys.org/ogs/ogs/gcc # Downloads the image to gcc_latest.sif
# OR: Pull the image docker://registry.opengeosys.org/ogs/ogs/gcc-gui image for compiling the Data Explorer
singularity shell gcc_latest.sif
[Now inside the container]
mkdir build; cd build
cmake ../ogs -DCMAKE_BUILD_TYPE=Release -DOGS_DISABLE_CCACHE=ON # OR set env var CCACHE_DIR
ninja
./bin/ogs
```

### Run OGS inside a Container (called from outside)

Once ogs executable is built it can be called from **outside** the container:

```bash
singularity exec gcc_latest.sif build/bin/ogs some/path/project.prj
```

This starts the container, mounts your home directory inside the container, passes the current working directory and runs the ogs executable (which is in your home directory which is mounted inside the container) with the passed project file. Everything works as expected and is transparent to the user. When ogs finishes the container stops and you returns to the host system.

## Container generator

You can download a prebuilt container from Docker Hub as shown above (e.g. `singularity pull docker://ogs6/gcc`). But we also provide a [container generator]({{< ref "container.md" >}}) to create a specific container for your needs.
