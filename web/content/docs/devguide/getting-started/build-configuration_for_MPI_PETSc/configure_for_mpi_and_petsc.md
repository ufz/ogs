+++
title = "Build configuration for MPI and PETSc"
date = "2022-05-11T09:48:39+02:00"
author = "Wenqing Wang"
weight = 5

[menu]
  [menu.devguide]
    parent = "getting-started"
+++

In OGS, the domain decomposition method (DDC) with MPI is used for parallel
 computing, and the [PETSc](https://petsc.org) package is employed to manage
 the partitioned vectors and matrices of the assembled global system of linear
 equations and to solve the assembled global system of linear equations as well.

To compiled OGS with PETSc, MPI c++ compiler wrappers like OpenMPI, MPICH (
has many derivatives, e.g intelMPI), and PETSc have to be installed as prerequisites.

## Install MPI

MPI c++ compiler wrappers are essentially installed on most of Linux clusters
 and supercomputers, and they are usually managed by the environment module of operating system.
What you have to do on that platforms is just to load a proper MPI wrapper module.
On Linux powered laptop, PC, or work station, you can install OpenMPI or MPICH via
 the package manager of the Linux distribution, e.g. RPM of RedHat, pacman of Arch.

## Install PETSc

If PETSc is not installed on your system, you can download the PETSc source code
 (<https://petsc.org/release/download/>), compile it
 and install by yourself according to the PETSc's installation instruction on
<https://petsc.org/release/install/>.

Here is an example of PETSc configuration:

```bash
./configure PETSC_ARCH=linux-fast  COPTFLAGS='-O3  --prefix=/home/me/opt/petsc --with-debugging=0 -march=native -mtune=native' CXXOPTFLAGS='-O3 -march=native -mtune=native' --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --download-fblaslapack --download-metis --download-parmetis --download-superlu_dist --download-scalapack --download-mumps  --download-hypre --with-c2html=0  --with-cxx-dialect=C++11 --with-cuda=0`
```

With that configuration, PETSc will be compiled without PETSc debugging model (`--with-debugging=0`),
 and installed to `opt` directory of the home directory. For Arch Linux, if you
 install PETSc package from AUR repository, please keep in mind that
 its option of `--with-debugging` is not set so far. This mean that PETSc is installed
 with PETSc debugging model by default.
 PETSc debugging model gives low performance. You can manually edit PKGBUILD of PETSc from AUR by
 adding the option of `--with-debugging=0` before installation.

Please note that the PETSc package is preferred to build in release model regardless of
 OGS being built with release or debug preset because that you may not like to debug the
 source code of PETSc and want to have a better computational performance.

PESTSc is recommended to use on Linux. On Windows, it only runs on
 UNIX emulator Cygwin. A detailed description about how to use PETSc on Windows
  is available on this PETSc site: <https://petsc.org/main/install/windows/>.

One the frontends of EVE cluster of UFZ, your can load specified MPI, PETSc and
 other modules by command:

```bash
source [path to ogs source code]/scripts/env/eve/petsc.sh
```

## Build configuration of OGS

For CMake configuration, the option of `OGS_USE_PETSC` must be set to true.
 It is `-DOGS_USE_PETSC=ON` for command line.

MPI C/C++ compiler has to be used for compilation. `OGS_USE_MPI` is set
on automatically once OGS_USE_PETSC is on.  On the platforms that the
software is not managed by environment module, the standard C/C++ compiler
is the default compiler. For such platforms, C/C++ compiler is better to
 be explicitly specified as MPI C/C++ compiler in CMake configuration in case
 the standard C/C++ compiler being used. If standard C/C++ compiler is used
when `OGS_USE_PETSC` is switched on, the compilation will fail with missed
library error.

For simplicity, you can first run CMake with MPI C/C++ compiler specified. For example:

```bash
cmake [path to source code] [-GNinja] -DOGS_USE_PETSC=ON -DCMAKE_C_COMPILER=/usr/bin/mpicc -DCMAKE_CXX_COMPILER=/usr/bin/mpic++
```

and then use `ccmake` or *cmake-gui* to configure further more CMake options.

The other basic build configurations are the same as that described in [Build configuration]({{<ref
"build-configuration">}}).

Once the CMake configuration is done, you can run `make` or `ninja` (if `-GNinja` is set)
 to compile the source code as that for the compilation for serial computing code.
