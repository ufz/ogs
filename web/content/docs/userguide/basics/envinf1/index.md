+++
date = "2018-02-27T11:00:13+01:00"
title = "Eve"
author = "Lars Bilke"
weight = 2
draft = true

aliases = [ "/docs/quickstart/basics/envinf1" ]

[menu]
  [menu.userguide]
    parent = "basics"
+++

## Introduction

Members of the Department Environmental Informatics of the Helmholtz Centre for Environmental Research - UFZ can use the `frontend1` and `frontend2`-machines which are tightly connected to the Eve cluster system.

## Select OGS version

You select an OGS version by loading a module:

```bash
module use /global/apps/modulefiles

# Examples:
module load ogs            # Loads latest (maybe unstable) ogs in standard config
module load ogs/head/petsc # Loads latest (maybe unstable) petsc config
module load ogs/6.0.9      # Loads stable version 6.0.9 in standard config, not released yet
```

You can select only one version at a time. Run `module purge` to unload all previously loaded modules.

See [Quickstart]({{< ref "introduction.md" >}}) for running instructions.
