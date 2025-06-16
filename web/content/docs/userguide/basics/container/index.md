+++
date = "2018-11-14T15:00:13+01:00"
title = "Running OGS in a container"
author = "Lars Bilke"
weight = 31
+++

This page describes how to **run** OGS with the help of a Linux container (for **users**) with [Apptainer](https://apptainer.org), formerly [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html).

## Prerequisites

* A running installation of Apptainer:
  * Available on Eve (UFZ),  all envinfs (UFZ), Juwels (SC Jülich), Taurus (TU Dresden).
  * See the Apptainers Admin guide for [install instructions](https://apptainer.org/docs/admin/latest/admin_quickstart.html#installation).

## Get a container image

### Recommended option on EVE: use prebuilt images

On the EVE cluster system (UFZ) you can use prebuilt images which can be used for [PETSc-based cluster jobs]({{< ref "parallel_computing_mpi.md#2a-use-a-container-to-launch-mpi-ogs" >}}).

### Option: Download a release image (preferred)

Simply download an image from the [releases]({{< ref "/releases" >}}) page.

### Option: Download image from the latest master-branch build

Simply download an image from the latest master-branch build:

<!-- vale off -->
* [ogs-serial.squashfs](https://vip.s3.ufz.de/ogs/public/container/ogs/master/ogs-serial.squashfs)
* [ogs-mkl.squashfs](https://vip.s3.ufz.de/ogs/public/container/ogs/master/ogs-mkl.squashfs)  (with and MKL Pardiso-support )
* [ogs-petsc.squashfs](https://vip.s3.ufz.de/ogs/public/container/ogs/master/ogs-petsc.squashfs) (with PETSC-support)
* [ogs-petsc-mkl.squashfs](https://vip.s3.ufz.de/ogs/public/container/ogs/master/ogs-petsc-mkl.squashfs) (with PETSC- and MKL Pardiso-support )
<!-- vale on -->

## Run OGS inside a Container (called from outside)

```bash
# Linux only:
apptainer exec ogs-6.x.x-serial.squashfs ogs some/path/project.prj
```

This starts the container, mounts your home directory inside the container, passes the current working directory and runs the `ogs` executable (in your home directory which is mounted inside the container) with the passed project file. When using
containers, everything stays transparent to the user. When OGS finishes the container stops and returns to the host system.

You can also specify the full executable path in the container:

```bash
# Works on macOS too:
apptainer exec ogs-6.x.x-serial.squashfs /usr/local/ogs/bin/ogs ...
```

Running a benchmark:

```bash
# Create output directories
mkdir -p _out _out_mpi
# Run serial benchmark
apptainer exec ogs-6.x.x-serial.squashfs ogs -o _out [ogs-sources]/Tests/Data/Mechanics/Linear/disc_with_hole.prj
# Run serial benchmark with output validation (via vtkdiff)
apptainer exec ogs-6.x.x-serial.squashfs ogs -o _out -r [ogs-sources]/Tests/Data/Mechanics/Linear [ogs-sources]/Tests/Data/Mechanics/Linear/disc_with_hole.prj
# Run parallel benchmark with MPI
mpirun -np 4 apptainer ogs ogs-6.x.x-openmpi-2.1.2.squashfs ogs -o _out_mpi [ogs-sources]/Tests/Data/Mechanics/Linear/disc_with_hole.prj
```

You can run other contained executables as well, e.g. `vtkdiff`:

```bash
apptainer exec ogs-6.x.x-serial.squashfs vtkdiff --help
```

You can interactively explore the container with `apptainer shell` (you can see that you are **in** the container because of the `Apptainer>` prefix of the shell):

```bash
# Shell into container
apptainer shell ogs-6.x.x-serial.squashfs
# List binaries in the container
Apptainer> ls $PATH
... ogs vtkdiff
# Exit the container and get back to your hosts shell
Apptainer> exit
```

## Custom Python environment for the container

For certain benchmarks or tutorials you may need additional Python packages. You can create a Python [virtual environment](https://virtualenv.pypa.io/en/latest/) inside the container (stored on your host) and install packages via `pip` (inside the container):

```bash
mkdir my-working-directory && cd my-working-directory
apptainer shell my-container.squashfs

# Now in the container
virtualenv .venv
source .venv/bin/activate
pip install numpy # or whatever you need
# run ogs
exit

# Now outside the container
# The virtualenv-directory .venv still persists
# If you want to run something in the container with exec, source the venv before:
apptainer exec my-container.squashfs bash -c 'source .venv/bin/activate && ogs ...'
```
