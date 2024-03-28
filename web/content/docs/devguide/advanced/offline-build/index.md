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

- The OGS source code. Just archive the full OGS source code directory (also containing the git repository in `.git`) and un-archive on the target system.
- The [CPM]({{< ref "cpm.md" >}}) source cache. It can be obtained via the [OGS package registry](https://gitlab.opengeosys.org/ogs/ogs/-/packages/) (see below).
- Optional: The external dependencies (for MFront, PETSc or LIS) source cache. It can be obtained via the [OGS package registry](https://gitlab.opengeosys.org/ogs/ogs/-/packages/)

## CPM

The CPM source cache may change over time. To get the required file id check the `package_file_id` field in `web/data/versions.json`, e.g. with:

```bash
$ jq -r '.cpm.package_file_id' web/data/versions.json
{{< dataFile "versions.cpm.package_file_id" >}} # <-- current version on master
```

Now simply download the file with:

```bash
curl https://gitlab.opengeosys.org/ogs/ogs/-/package_files/[insert ID here]/download --output cpm.tar.gz
```

Un-archive the CPM cache into a directory. Configure OGS as usual but point to the extracted CPM cache:

```bash
cmake --preset release -DCPM_SOURCE_CACHE=./path/to/cpm
```

There will be some CMake warnings from CPM regarding missing git repositories in the cache. You can ignore them.

## External dependencies

The external dependencies source cache may change over time. To get the required version check the `cache_hash` field in `web/data/versions.json`, e.g. with:

```bash
$ jq -r '.ext.cache_hash' web/data/versions.json
{{< dataFile "versions.ext.cache_hash" >}} # <-- current version on master
```

On the [external dependencies package page](https://gitlab.opengeosys.org/ogs/ogs/-/packages) download the file `ext.tar.gz` of the package `external-dependencies` with the specified version.

Extract the external dependencies cache into a directory. Configure OGS as usual but point to the extracted external dependencies cache:

```bash
cmake --preset release -DOGS_EXTERNAL_DEPENDENCIES_CACHE=./path/to/ext
```

## Example for OpenGeoSys 6.5.1

On a machine with internet access download the tarballs:

```bash
wget https://gitlab.opengeosys.org/ogs/ogs/-/archive/6.5.1/ogs-6.5.1.tar.gz
wget -O cpm.tar.gz https://gitlab.opengeosys.org/ogs/ogs/-/package_files/1048/download
wget -O ext.tar.gz https://gitlab.opengeosys.org/ogs/ogs/-/package_files/1245/download
```

Copy those file onto the machine where you want to build OGS. Then build OGS offline:

```bash
tar xf ogs-6.5.1.tar.gz
tar xf cpm.tar.gz
tar xf ext.tar.gz

cd ogs-6.5.1
OGS_VERSION=6.5.1 cmake --preset release -DCPM_SOURCE_CACHE=../cpm -DOGS_EXTERNAL_DEPENDENCIES_CACHE=../ext
cmake --build --preset release
```
