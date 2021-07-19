+++
date = "2018-02-26T11:00:13+01:00"
title = "Conan package manager"
author = "Lars Bilke"
weight = 1042

aliases = ["/docs/devguide/advanced/conan-package-manager"]

[menu]
  [menu.devguide]
    parent = "packages"
+++

<div class='note'>

### <i class="far fa-exclamation-triangle"></i> Conan {{< dataFile "versions.minimum_version.conan" >}} required

A recent version ({{< dataFile "versions.minimum_version.conan" >}}) of Conan is required! Please update Conan by running `pip install --upgrade conan` or by downloading the Windows installer.
</div>

The [Conan package manager](https://www.conan.io) helps to install all required libraries in a convenient way on every platform. See [Setup pre-requisites](../../getting-started/prerequisites) for installation instructions. If the Conan executable is found Conan is used for third-party library handling. Set the CMake option `OGS_USE_CONAN=OFF` to disable Conan.

## Advanced usage

### Build packages locally

Per default when Conan is enabled it will try to fetch prebuilt binaries from the [OGS Conan repository](https://ogs.jfrog.io/ogs/conan/) at <https://ogs.jfrog.io/ogs/api/conan/conan>. With the CMake option `OGS_CONAN_BUILD` you define what gets build locally. This option can be set to:

- `missing` - Default, only builds packages which are not available as a prebuilt binary for the current configuration
- `all` - Builds all packages locally
- `never` - Builds no package locally
- `[a list of libraries to build]`, e.g. `"petsc;tfel"`. For names see [ConanSetup.cmake](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/scripts/cmake/ConanSetup.cmake). Make sure to set this back to 'never' after the libs have been built. Otherwise it would rebuild the libs on the next CMake run.

### Conan environment

Conan creates a [runtime environment](https://docs.conan.io/en/latest/mastering/virtualenv.html#virtualrunenv-generator) automatically. When activated this sets environment variables such as `PATH` or `LD_LIBRARY_PATH` to point to the used Conan packages. Can be activated in your build directory after CMake ran:

```bash
source activate_run.sh
...
source deactivate_run.sh
```

## Further information

- [Conan Blog](https://blog.conan.io)
- [Conan Documentation](https://docs.conan.io/en/latest/)
- [Bincrafters Blog](https://bincrafters.github.io/tag/conan/), a blog from a group of active Conan package developers
- [Bincrafters Documentation](https://bincrafters.readthedocs.io/en/latest/introduction_to_conan.html)
