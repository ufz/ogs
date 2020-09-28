+++
date = "2020-09-28"
title = "Run"
author = "Lars Bilke"
weight = 1045

[menu]
  [menu.devguide]
    parent = "troubleshooting"
    identifier = "run-troubleshooting"
+++

This page describes errors you get at runtime of OGS, e.g. when executing the `ogs`-executable or some utilities.

## Error message: `error while loading shared libraries` / `XX.dll could not be found or libXX.so could not be found.`

### Linux / macOS

Typical error message:

```
error while loading shared libraries: libXX.so: cannot open shared object file: No such file or directory
```

A shared library which was linked to OGS could not be found during runtime. The runtime search paths are determined by the system configuration but you add paths with the environment variable `LD_LIBRARY_PATH` (macOS: `DYLD_LIBRARY_PATH`). So you know where the missing library is located you can adapt the environment variable:

```bash
export LD_LIBRARY_PATH=/path/to/missing/lib:$LD_LIBRARY_PATH
./ogs
```
More information can be found here:

- https://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html
- https://amir.rachum.com/blog/2016/09/17/shared-libraries/#runtime-search-path

If you use Conan you want to activate its [runtime environment]({{< relref "conan-package-manager.md#conan-environment" >}}).

One can see with the `ldd` (macOS: `otool -L`) tool which dynamic libraries will be loaded at runtime.

### Windows

Typical error message:

```
XX.dll could not be found.
```

Similar to Linux but on Windows you need to adapt the `PATH` environment variable.

One can use the [Dependencies](https://github.com/lucasg/Dependencies)-tool to see which DLLs will be loaded at runtime.
