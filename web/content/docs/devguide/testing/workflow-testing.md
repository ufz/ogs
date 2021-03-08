+++
date = "2021-03-05T10:49"
title = "Workflow testing"
author = "Lars Bilke"
weight = 1024

[menu]
  [menu.devguide]
    parent = "testing"
+++

## Introdution

We use the workflow manager [Snakemake](https://snakemake.readthedocs.io) to test workflows which consist of the execution of several steps which are based on each other.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) >= {{< dataFile "versions.minimum_version.snakemake" >}}
    - If you use [Poetry]({{< ref "python-env.md#poetry" >}}) then `snakemake` is installed in your virtual environment in your build-directory automatically. You can then call it via `poetry run snakemake ...`.
- On Windows only:
    - The `tee`-utility in the `PATH` (can be installed from https://sourceforge.net/projects/unxutils)

## Examples

- [ExtractBoundary.smk](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Applications/Utils/ExtractBoundary.smk)
- [VoxelGridFromLayers.smk](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Applications/Utils/VoxelGridFromLayers.smk)

These example workflows [are added to ctest](https://gitlab.opengeosys.org/ogs/ogs/-/blob/540d0b454c9e3805a81f7c4a1b6ee7565be6845c/Applications/Utils/Tests.cmake#L302-315) as well:

```cmake
if(SNAKEMAKE AND NOT OGS_USE_MPI)
    add_test(NAME snakemake_ExtractBoundary
        COMMAND ${SNAKEMAKE} -j 1
            --configfile ${PROJECT_BINARY_DIR}/buildinfo.yaml
            -s ${CMAKE_CURRENT_SOURCE_DIR}/ExtractBoundary.smk
    )
    add_test(NAME snakemake_VoxelGridFromLayers
    # ...
    )
    add_dependencies(ctest ExtractBoundary Layers2Grid AddFaultToVoxelGrid)
endif()
```

## Modularization

We started on implementing modular rule definitions and tool wrapper in [`scripts/snakemake`](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/scripts/snakemake).

## Links

- [Snakemake Documentation](https://snakemake.readthedocs.io)
- [Short tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/short.html)
- [Software Carpentry Workshop](https://carpentries-incubator.github.io/workflows-snakemake/index.html)
- [HPC Carpentry Workshop with Snakemake](http://www.hpc-carpentry.org/hpc-python/11-snakemake-intro/index.html)
