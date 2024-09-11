+++
title = "Build configuration for MPI and PETSc"
date = "2022-05-11T09:48:39+02:00"
author = "Wenqing Wang, Feliks Kiszkurno, Lars Bilke"
weight = 5

[menu]
  [menu.devguide]
    parent = "getting-started"
+++

In OGS, the domain decomposition method (DDC) with MPI is used for parallel computing, and the [PETSc](https://petsc.org) package is employed to manage the partitioned vectors and matrices of the assembled global system of linear equations and to solve the assembled global system of linear equations as well.

To compile OGS with PETSc, MPI c++ compiler wrappers like OpenMPI, MPICH (
has many derivatives, e.g intelMPI) has to be installed as prerequisite.

<div class='win'>

<div class='note'>

PETSc is not supported on Windows system.
As an alternative you can follow the instructions in Section [Install PETSc manually](#install-petsc-manually) and run PETSc using Cygwin or use Windows Subsystem For Linux (WSL).
The latter option is recommended - using WSL.
The manual of setting up OpenGeoSys in WSL can be found in our [Windows Subsystem For Linux]({{< ref "wsl" >}}) guide.
After setting up WSL, please follow the Linux tab in this guide.

</div>

</div>

## Set up prerequisites

Before continuing with this guide, please follow all steps from the "Developer guide" articles: [Set Up Prerequisites]({{< ref "prerequisites" >}}) and [Get the source code]({{< ref "get-the-source-code" >}}).

### Install `pkg-config`

Is used for finding the PETSc installation. Can be installed via e.g. `apt` or `brew`:

```bash
apt-get install pkg-config # Linux
brew install pkg-config    # macOS
```

### Install MPI

<div class='win'>

Please setup OpenGeoSys environment in WSL and proceed with instructions from Linux tab (recommended), or follow the Section [Install PETSc manually](#install-petsc-manually) from Windows tab (not recommended).

</div>

<div class='linux'>

MPI c++ compiler wrappers are essentially installed on most of Linux clusters and supercomputers, and they are usually managed by the environment module of operating system.
What you have to do on that platforms is just to load a proper MPI wrapper module.
On Linux powered laptop, PC, or work station, you can install OpenMPI or MPICH via the package manager of the Linux distribution, e.g. apt in Ubuntu, RPM in RedHat, pacman in ArchLinux.

On Ubuntu it is sufficient to install the package `libopenmpi-dev`.

Installation can be done using following command:

```bash
sudo apt install libopenmpi-dev
```

On other distributions, the exact names of the packages can slightly differ.
If this is the case, you can search for the correct names using the search function of your distributions package manager.

</div>

<div class='mac'>

On MacOS machines, you can install OpenMPI using ```brew``` package manager:

```bash
brew install open-mpi
```

</div>

## Build configuration of OGS

### Automatic way

The easiest way of obtaining a working copy of OpenGeoSys with PETSc and MPI, is to use the ```release-petsc``` preset.
It will automatically download and compile a version of PETSc compatible with OpenGeoSys (if not system installed PETSc is found).
This could be done by running following commands inside of the folder containing OpenGeoSys source code:

```bash
cmake --preset release-petsc
```

and

```bash
cmake --build --preset release-petsc
```

The ```release-petsc``` preset is based on the ```release``` preset.
It can be customized by passing additional variable with ```-D``` flags to `cmake` command.

<div class='win'>
</div>
<div class='linux'>
</div>
<div class='mac'>

<div class='note'>

### <i class="far fa-exclamation-triangle"></i> Additional remarks for Homebrew users

If you use the `homebrew`-package manager please be aware that if you have installed the following packages:

- `hdf5`
- `vtk`

Then you need to add these CMake-variables to your configuration:

- `-DOGS_BUILD_HDF5=ON`
- `-DOGS_BUILD_VTK=ON`

Otherwise the `brew`-installed libraries are picked up which are built for the serial case and are therefore incompatible with the PETSc configuration of OGS.

When building VTK locally you need to also install [git-lfs](https://git-lfs.com) via:

```bash
brew install git-lfs
```

</div>

</div>

If you want to install PETSc yourself, please see Section [Install PETSc](#install-petsc-manually).

### Manual way

If settings of the preset ```release-petsc``` doesn't suit your use case and you would like to enable PETSc in other preset, please follow the instructions in this section.

For CMake configuration, the option of `OGS_USE_PETSC` must be set to true.

A MPI C/C++ compiler has to be used for compilation.
You can specify it via the environment variables `CC` and `CXX`, e.g.:

```bash
CC=mpicc CXX=mpic++ cmake --preset my-preset [...]
```

The other basic build configurations are the same as that described in [Build configuration]({{<ref"build-configuration">}}).

Once the CMake configuration is done, you can run `ninja` (or `make`) to compile the source code, same as compilation for the non-PETSc build.

## Install PETSc manually

The instructions from this section apply only if you want to install or compile PETSc yourself, instead of letting OpenGeoSys do it automatically.

<div class='linux'>

If PETSc is not installed on your system, you can download the PETSc source code
 (<https://petsc.org/release/download/>), compile it and install by yourself according to the PETSc's installation instruction on <https://petsc.org/release/install/>.

Here is an example of PETSc configuration:

```bash
./configure PETSC_ARCH=linux-fast  COPTFLAGS='-O3  --prefix=/home/me/opt/petsc --with-debugging=0 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --download-fblaslapack --download-metis --download-parmetis --download-superlu_dist --download-scalapack --download-mumps  --download-hypre --with-c2html=0  --with-cxx-dialect=C++11 --with-cuda=0`
```

With that configuration, PETSc will be compiled with optimizations and without debugging info (`--with-debugging=0`), and installed to user's `opt` directory.
For Arch Linux, if you install PETSc package from AUR repository, please keep in mind that its option of `--with-debugging` is not set so far.
You can manually edit PKGBUILD of PETSc from AUR by adding the option of `--with-debugging=0` before installation.

Please note, that the PETSc package is usually build in release mode independent of the OGS' release or debug build.

</div>

<div class='win'>

PETSc is recommended to use on Linux.
On Windows, it only runs on UNIX emulator Cygwin.
A detailed description about how to use PETSc on Windows is available on this PETSc site: <https://petsc.org/release/install/windows/>.

</div>

## Using PETSc on EVE cluster (UFZ)

On the frontends of the UFZ EVE cluster the MPI, PETSc, and other required modules can be loaded with the following command:

```bash
source [path to ogs source code]/scripts/env/eve/petsc.sh
```
