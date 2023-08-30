+++
date = "2018-02-26T11:00:13+01:00"
title = "Working on Eve"
author = "Lars Bilke"
weight = 1053

aliases = ["/docs/devguide/advanced/working-on-eve"]

[menu]
  [menu.devguide]
    parent = "environments"
+++

## Introduction

Members of the Department Environmental Informatics of the Helmholtz Centre for Environmental Research - UFZ can use the `frontend1`- and `frontend2`-machines from the Eve cluster system.

## Run OGS within a container

Eve has the [Apptainer/Singularity](https://apptainer.org) container runtime installed so you can simply download a prebuilt container and run OGS inside it. See the [user guide]({{< ref "container" >}}) for more information.

## Build OGS-6

Load required modules by sourcing the environment script:

```bash
source scripts/env/eve/cli.sh
```

Then do the [build configuration]({{< ref "build-configuration.md" >}}) and [build]({{< ref "/docs/devguide/getting-started/build.md" >}}) the project.
