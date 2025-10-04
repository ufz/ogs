+++
date = "2018-02-26T11:00:13+01:00"
title = "Redistributable builds"
author = "Lars Bilke"
weight = 1052

aliases = ["/docs/devguide/advanced/redistributable-builds"]

[menu]
  [menu.devguide]
    parent = "environments"
+++

## Introduction

If a binary runs on a different machine depends on a lot of factors. The following sections provide some tips to maximize portability. In general Windows binaries are more portable than *nix binaries.

## CPU Architecture Optimization

<div class='win'>

### Default compiler flags

- `/favor:blend`: Optimizes for AMD **and** Intel architectures, see [MSDN docs](https://msdn.microsoft.com/en-us/library/ms173505.aspx)
- `/arch:sse2` on [32-bit](https://msdn.microsoft.com/en-us//library/7t5yh4fd.aspx), `/arch` undefined on [64-bit](https://msdn.microsoft.com/en-us/library/jj620901.aspx)

### Optimization CMake options

- `OGS_CPU_ARCHITECTURE`: possible values `AMD64` or `INTEL64` (sets `/favor:{AMD64 | INTEL64}` flag)

</div>

<div class='linux'>

### Default compiler flags

- `-march=native`: Optimizes for current CPU

### CMake options

#### For redistribution

- `OGS_CPU_ARCHITECTURE`: set to `generic` for good balance between optimization and portability; set to `core2` for maximum portability, more info on [GCC docs](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html)

#### For optimization

- `OGS_CPU_ARCHITECTURE`: tot to `native` for best optimization for your current CPU, possible values are [listed here](https://stackoverflow.com/a/25095818/80480), more info on [GCC docs](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html)

</div>

<div class='mac'>
See Linux-tab!
</div>

## Package the build

Use the `package`-target which tries to gather all dependencies and fixes up shared library paths:

```bash
make package
```

This creates a zip- or tar-archive which should be redistributable.
