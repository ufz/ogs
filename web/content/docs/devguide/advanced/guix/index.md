+++
date = "2023-07-20"
title = "Build with GNU Guix"
author = "Lars Bilke"
weight = 1069

[menu]
  [menu.devguide]
    parent = "advanced"
+++

<div class='note'>

### <i class="far fa-info-circle"></i> Warning

This page and ogs builds with GNU Guix are currently work-in-progress!

</div>

[GNU Guix](https://guix.gnu.org) is a distribution of the GNU operating system but also a package manager. You can use it to create bit-by-bit reproducible (Linux) binaries of `ogs`.

The package definitions for OGS are defined in [this repo](https://gitlab.opengeosys.org/ogs/inf/guix-ogs) which can be used as a Guix channel. `guix` is currently installed on `envinf3`.

## Building

### From local source tree

```bash
# builds ogs serial config
guix time-machine -C scripts/guix/channels.scm -- \
  build ogs-serial --with-source=ogs@6.4.4-testing=$PWD
# SteadyStateDiffuion process only:
guix time-machine -C scripts/guix/channels.scm -- \
  build ogs-ssd --with-source=ogs-ssd@6.4.4-testing=$PWD
# ogs petsc config
guix time-machine -C scripts/guix/channels.scm -- \
  build ogs-petsc --with-source=ogs-petsc@6.4.4-testing=$PWD
```

The version string starting after `@` can be set individually. To get the dependency tree:

```bash
guix time-machine -C scripts/guix/channels.scm -- graph ogs-serial | dot -Tpdf > dag.pdf
```

### From web URL

```bash
wget https://gitlab.opengeosys.org/ogs/ogs/-/raw/master/scripts/guix/channels.scm
guix time-machine -C ./channels.scm -- build ogs \
  --with-source=ogs@6.4.4-testing=https://gitlab.opengeosys.org/ogs/ogs/-/archive/master/ogs-master.tar.bz2
```

## Developing

```bash
guix time-machine -C scripts/guix/channels.scm -- \
  shell --container --nesting --network --development ogs-serial \ # OR ogs-petsc
  openssl nss-certs coreutils bash git
# Now in guix shell with all dependencies installed:
export CMAKE_PRESET_BUILD_DIR_PREFIX=guix/ # presets then create e.g. ../build/guix/release
cmake --preset release -DOGS_BUILD_PROCESSES=SteadyStateDiffusion
cmake --build --preset release
```

As a shortcut you can use this script:

```bash
./scripts/guix/ogs-env.sh ogs-serial # or `ogs-petsc` for parallel
```

You can also play with e.g. different versions of dependencies, here changing OpenMPI to 4.1.6:

```bash
guix time-machine -C scripts/guix/channels.scm -- \
  shell --container --nesting --network --development ogs-petsc \
  openssl nss-certs coreutils bash git \
  --with-source=openmpi@4.1.6=https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.6.tar.bz2
```

## Links

- [`guix build`](https://guix.gnu.org/manual/en/html_node/Invoking-guix-build.html)-reference
- [`guix shell`](https://guix.gnu.org/manual/en/html_node/Invoking-guix-shell.html)-reference
- [Blog: The ultimate guide to software development with Guix](https://guix.gnu.org/en/blog/2023/from-development-environments-to-continuous-integrationthe-ultimate-guide-to-software-development-with-guix/)
- [Futurile blog](https://www.futurile.net)
- [Peter's Blog](https://peterloleungyau.github.io/post/)
