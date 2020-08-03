+++
date = "2018-09-21T11:00:13+01:00"
title = "Singularity"
author = "Lars Bilke"
weight = 1036

[menu]
  [menu.devguide]
    parent = "advanced"
+++

## Introduction

[Singularity](https://www.sylabs.io) is a Linux container runtime similar to Docker. Key advantages over Docker are

- Container don't run with root privileges
- You are the same user with the same privileges inside the container as on the host
- Container can run on HPC systems and seamlessly integrate with resource managers and MPI
- Container can leverage NVidia GPUs

Singularity per default mounts your home directory and also passes your current working directory when starting a container. Therefore it is easy to use it for development.

### Prerequisites

- Linux and [Mac](https://sylabs.io/singularity-desktop-macos/) only!
  - Note: You can run Singularity on Windows inside a virtual machine using [WSL 2](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
- [Install Git]({{< ref "prerequisites" >}}/#step-install-git)
- [Install Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps)

### Build OGS inside a container

```bash
[git clone ogs]
singularity pull singularity pull docker://registry.opengeosys.org/ogs/ogs/ogs/gcc:master # Downloads the image to gcc.sif
# OR: Pull the singularity pull docker://registry.opengeosys.org/ogs/ogs/ogs/gcc-gui:master image for compiling the Data Explorer
singularity shell gcc_master.sif
[Now inside the container]
mkdir build; cd build
cmake ../ogs -DCMAKE_BUILD_TYPE=Release
make
./bin/ogs
```

### Run OGS inside a Container (called from outside)

Once ogs executable is built it can be called from **outside** the container:

```bash
singularity exec gcc_master.sif build/bin/ogs some/path/project.prj
```

This starts the container, mounts your home directory inside the container, passes the current working directory and runs the ogs executable (which is in your home directory which is mounted inside the container) with the passed project file. Everything works as expected and is transparent to the user. When ogs finishes the container stops and you returns to the host system.

## Container generator

You can download a prebuilt container from Docker Hub as shown above (e.g. `singularity pull docker://ogs6/gcc`). But we also provide a [container generator]({{< ref "container.md" >}}) to create a specific container for your needs.
