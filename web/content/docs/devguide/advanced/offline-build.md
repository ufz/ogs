+++
date = "2022-02-09T11:00:13+01:00"
title = "Offline build"
author = "Lars Bilke"
weight = 1068

[menu]
  [menu.devguide]
    parent = "advanced"
+++

OGS can be built on systems without internet connection when the following files can be made available on the system:

- The ogs source code. Just archive the full ogs source code directory (also containing the git repository in `.git`) and unarchive on the target system.
- The [CPM]({{< ref "cpm.md" >}}) source cache. It can be obtained via the [OGS package registry](https://gitlab.opengeosys.org/ogs/ogs/-/packages/) (see below).
- Optional: The external dependencies (for MFront, PETSc or LIS) source cache. It can be obtained via the [OGS package registry](https://gitlab.opengeosys.org/ogs/ogs/-/packages/)

## CPM

The cpm source cache may change over time. To get the required version check the `cache_hash` field in `web/data/versions.json`, e.g. with:

```bash
$ jq -r '.cpm.cache_hash' web/data/versions.json
{{< dataFile "versions.cpm.cache_hash" >}} # <-- current version on master
```

On the [cpm package page](https://gitlab.opengeosys.org/ogs/ogs/-/packages/1) download the file `cpm.tar.gz` with the specified version.

Unarchive the cpm cache into a directory. Configure OGS as usual but point to the extracted cpm cache:

```bash
cmake -S . --preset release -DCPM_SOURCE_DIR=./path/to/cpm
```

There will be some CMake warnings from CPM regarding missing git repositories in the cache. You can ignore them.

## External dependencies

The external dependencies source cache may change over time. To get the required version check the `cache_hash` field in `web/data/versions.json`, e.g. with:

```bash
$ jq -r '.ext.cache_hash' web/data/versions.json
{{< dataFile "versions.ext.cache_hash" >}} # <-- current version on master
```

On the [external dependencies package page](https://gitlab.opengeosys.org/ogs/ogs/-/packages/14) download the file `ext.tar.gz` with the specified version.

Unarchive the external dependencies cache into a directory. Configure OGS as usual but point to the extracted external dependencies cache:

```bash
cmake -S . --preset release -DOGS_EXTERNAL_DEPENDENCIES_CACHE=./path/to/ext
```
