+++
date = "2021-02-11T10:46"
title = "CMake dependency management"
author = "Lars Bilke"
weight = 1032

aliases = ["/docs/devguide/advanced/cpm"]

[menu]
  [menu.devguide]
    parent = "packages"
+++

We employ [CPM](https://github.com/cpm-cmake/CPM.cmake#options), a CMake dependency management solution, to integrate third-party dependencies.

[Dependencies](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/scripts/cmake/Dependencies.cmake) are downloaded at configure-time into the build directory. Dependencies get build in the build directory inside the `_deps`-subdirectory.

To speed things up [you can specify a cache directory](https://github.com/cpm-cmake/CPM.cmake#cpm_source_cache) for the downloads by setting `CPM_SOURCE_CACHE` either as an environment or CMake variable:

```bash
# Can be placed in .bashrc or .bash_profile:
export CPM_SOURCE_CACHE=$HOME/.cache/CPM
# OR
cmake ... -DCPM_SOURCE_CACHE=$HOME/.cache/CPM
```

Some dependencies (those which are added with `CPMFindPackage()`) are first searched to be locally installed on the system with a fallback to CPM if not found. You can disable the search for local packages with [`CPM_DOWNLOAD_ALL`](https://github.com/cpm-cmake/CPM.cmake#cpm_download_all).

## Further information

- [CPM on GitHub](https://github.com/cpm-cmake/CPM.cmake#options)
- [Blog post: An Awesome Dependency Manager for C++ with CMake](https://medium.com/swlh/cpm-an-awesome-dependency-manager-for-c-with-cmake-3c53f4376766)
- [Blog post: CMake and the Future of C++ Package Management](https://ibob.github.io/blog/2020/01/13/cmake-package-management/)
