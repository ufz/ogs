+++
date = "2018-02-26T11:00:13+01:00"
title = "Working on Eve"
author = "Lars Bilke"
weight = 1040

aliases = ["/docs/devguide/advanced/working-on-eve"]

[menu]
  [menu.devguide]
    parent = "environments"
+++

## Introduction

Members of the Department Environmental Informatics of the Helmholtz Centre for Environmental Research - UFZ can use the `frontend1`- and `frontend2`-machines from the Eve cluster system.

<div class='note'>
You need to opt-in to use the [EasyBuild modules](https://wiki.ufz.de/eve/index.php/EasyBuild)! To do so follow these steps:

- Login to `frontend1` or `frontend2`
- Run `touch ~/.easybuild-yes`
- `exit` and re-login

Also on eve the installed git module (2.10) is over 2 years old and we recommend using a newer module:

```bash
ml use /global/apps/modulefiles
ml git/2.20.1
```

</div>

## Run OGS within a container

Eve has [Singularity](https://www.sylabs.io/singularity) container runtime installed so you can simply download a prebuilt container and run OGS inside it. See the [user guide]({{< ref "container" >}}) for more infos.

## Build OGS-6

Load required modules by sourcing the environment script:

```bash
source scripts/env/eve/cli.sh
```

Then do the [build configuration]({{< ref "build-configuration.md" >}}) and [build]({{< ref "/docs/devguide/getting-started/build.md" >}}) the project.

## Optional steps

### Install and use Conan

Follow the instructions on using Python's `virtualenv` [on the Eve-wiki](https://wiki.ufz.de/eve/index.php/Python#virtualenv) for setting up a local Python environment. Then you can `pip install conan` and use Conan.

### Additional Features

Generating Doxygen documentation:

```bash
module load Doxygen/1.8.14
```

You can build with Ninja:

```bash
module load ninja/1.9.0
```

## Build OGS-5

```bash
module load CMake/3.11.4
module load GCC/6.4.0-2.28
```
