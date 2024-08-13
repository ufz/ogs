+++
date = "2018-11-14T15:00:13+01:00"
title = "Running OGS in a container"
author = "Lars Bilke"
weight = 3
+++

<div class='note'>

### Important note

<!-- TODO: Consider to place this section somewhere outside of basics, e.g., as an individual section as the context may be
considered already pretty advanced. -->
This page describes how to **run** OGS with the help of a Linux container (for **users**).
To **build** OGS with the help of a container go to the [Developer Guide]({{< ref "singularity.md" >}}) (for **developers**).

</div>

## With Singularity
<!-- TODO: Update this section regarding Apptainer -->

### Prerequisites

* A running installation of Singularity 3.0 or higher:
  * Available on Eve (UFZ), `envinf1` / `envinf2` (UFZ), Juwels (SC Jülich), Taurus (TU Dresden).
  * See the developer guide for [install instructions]({{< ref "singularity.md#prerequisites" >}}).

### Get a container image

#### Recommended option on EVE: use prebuilt images

On the EVE cluster system (UFZ) you can use prebuilt images which can be used for [PETSc-based cluster jobs]({{< ref "parallel_computing_mpi.md#2a-use-a-container-to-launch-mpi-ogs" >}}).

#### Option: Download a release image (preferred)

Simply download an image from the [releases]({{< ref "/releases" >}}) page.

#### Option: Download image from the latest master-branch build

Simply download an image from the latest master-branch build:

<!-- vale off -->
* [ogs-serial.squashfs](https://minio.ufz.de/ogs/public/container/ogs/master/ogs-serial.squashfs)
* [ogs-mkl.squashfs](https://minio.ufz.de/ogs/public/container/ogs/master/ogs-mkl.squashfs)  (with and MKL Pardiso-support )
* [ogs-petsc.squashfs](https://minio.ufz.de/ogs/public/container/ogs/master/ogs-petsc.squashfs) (with PETSC-support)
* [ogs-petsc-mkl.squashfs](https://minio.ufz.de/ogs/public/container/ogs/master/ogs-petsc-mkl.squashfs) (with PETSC- and MKL Pardiso-support )
<!-- vale on -->

Please note that these images have `.squashfs` as the file ending. Usage is the same as with `.sif`-files.

### Run OGS inside a Container (called from outside)

```bash
# Linux only:
singularity exec ogs-6.x.x-serial.sif ogs some/path/project.prj
```

This starts the container, mounts your home directory inside the container, passes the current working directory and runs the `ogs` executable (in your home directory which is mounted inside the container) with the passed project file. When using
containers, everything stays transparent to the user. When OGS finishes the container stops and returns to the host system.

You can also specify the full executable path in the container:

```bash
# Works on macOS too:
singularity exec ogs-6.x.x-serial.sif /usr/local/ogs/bin/ogs ...
```

Running a benchmark:

```bash
# Create output directories
mkdir -p _out _out_mpi
# Run serial benchmark
singularity exec ogs-6.x.x-serial.sif ogs -o _out [ogs-sources]/Tests/Data/Mechanics/Linear/disc_with_hole.prj
# Run serial benchmark with output validation (via vtkdiff)
singularity exec ogs-6.x.x-serial.sif ogs -o _out -r [ogs-sources]/Tests/Data/Mechanics/Linear [ogs-sources]/Tests/Data/Mechanics/Linear/disc_with_hole.prj
# Run parallel benchmark with MPI
mpirun -np 4 singularity ogs ogs-6.x.x-openmpi-2.1.2.sif ogs -o _out_mpi [ogs-sources]/Tests/Data/Mechanics/Linear/disc_with_hole.prj
```

You can run other contained executables as well, e.g. `vtkdiff`:

```bash
singularity exec ogs-6.x.x-serial.sif vtkdiff --help
```

You can interactively explore the container with `singularity shell` (you can see that you are **in** the container because of the `Singularity [container image file]:...>` prefix of the shell):

```bash
# Shell into container
singularity shell ogs-6.x.x-serial.sif
# List files in the container
Singularity ogs-6.x.x-serial.sif:...> ls /scif/apps/ogs/bin
... ogs tetgen vtkdiff
# Exit the container and get back to your hosts shell
Singularity ogs-6.x.x-serial.sif:...> exit
```

### Custom Python environment for the container

For certain benchmarks or tutorials you may need additional Python packages. You can create a Python [virtual environment](https://virtualenv.pypa.io/en/latest/) inside the container (stored on your host) and install packages via `pip` (inside the container):

```bash
mkdir my-working-directory && cd my-working-directory
singularity shell my-container.sif

# Now in the container
virtualenv .venv
source .venv/bin/activate
pip install numpy # or whatever you need
# run ogs
exit

# Now outside the container
# The virtualenv-directory .venv still persists
# If you want to run something in the container with exec, source the venv before:
singularity exec my-container.sif bash -c 'source .venv/bin/activate && ogs ...'
```

### Run the DataExplorer inside a Container

* Get a Singularity container with the DataExplorer (has `-gui` in its name)
* `singularity exec ogs-xxx-gui-xxx.sif DataExplorer`

You may use this container on e.g. `envinf1` with X11 forwarding (`ssh -XY envinf1`).

----

## With Docker

Although Singularity is the preferred container runtime you can use [Docker](https://www.docker.com) too.

### Prerequisites

* A running [installation of Docker](https://docs.docker.com/get-docker/)

### Run OGS inside a Docker container

* Get the container: `docker pull registry.opengeosys.org/ogs/ogs/ogs-serial`
* Start interactive container session: `docker run --rm -it registry.opengeosys.org/ogs/ogs/ogs-serial`
* Run `ogs`: `/usr/local/ogs/bin/ogs --version`
* Exit the container: `exit`

You will notice that the interactive session in your container is isolated from your host, i.e. you do not have access to files on your host. You need to explicitly [mount](https://docs.docker.com/storage/bind-mounts/) them on `docker run`:

```bash
mkdir ~/ogs_out
docker run --rm -it -v $HOME/code/ogs6/ogs/Tests/Data:/tmp/data:ro -v $HOME/ogs_out:/tmp/out registry.opengeosys.org/ogs/ogs/ogs-serial
/usr/local/ogs/bin/ogs -o /tmp/out /tmp/data/Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e4.prj
exit
ls ~/ogs_out
# [shows ogs generated output files]
```
