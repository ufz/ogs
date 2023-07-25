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

[GNU Guix](https://guix.gnu.org) is a distribution of the GNU operating system but also a package manager. You can use it to create bit-by-bit reproducible (Linux) binaries of `ogs`. The package definitions are defined in [`guix.scm`](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/guix.scm). `guix` is currently installed on `envinf3`.

## Building

```bash
guix build -f guix.scm                     # builds ogs serial config
guix build -L $PWD/.guix/modules ogs-ssd   # SteadyStateDiffuion process only
guix build -L $PWD/.guix/modules ogs-petsc # ogs petsc config
```

## Developing

```bash
guix shell
# Now in guix shell with all dependencies installed:
export CMAKE_PRESET_BUILD_DIR_PREFIX=guix/ # presets then create e.g. ../build/guix/release
cmake --preset release -DOGS_BUILD_PROCESSES=SteadyStateDiffusion
cmake --build --preset release
```

## Links

- [`guix build`](https://guix.gnu.org/manual/en/html_node/Invoking-guix-build.html)-reference
- [`guix shell`](https://guix.gnu.org/manual/en/html_node/Invoking-guix-shell.html)-reference
- [`guix install`](https://guix.gnu.org/manual/en/html_node/Invoking-guix-install.html)-reference
- [Blog: The ultimate guide to software development with Guix](https://guix.gnu.org/en/blog/2023/from-development-environments-to-continuous-integrationthe-ultimate-guide-to-software-development-with-guix/)
