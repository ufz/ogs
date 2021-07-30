+++
date = "2018-02-26T11:00:13+01:00"
title = "Using a compiler cache"
author = "Lars Bilke"
weight = 1062

aliases = ["/docs/devguide/advanced/using-ccache"]

[menu]
  [menu.devguide]
    parent = "advanced"
+++

## Introduction

A compiler cache speeds up compilation times by caching object files and reusing them on subsequent builds. This even works for complete rebuilds (i.e. deleting the full build-directory). The compiler cache [ccache](https://ccache.samba.org) is automatically used when it is found by CMake.

<div class='linux'>

Install it with your package manager, e.g.:

```bash
sudo apt install ccache
```

## Usage on eve

Just load the module:

```bash
module load /global/apps/modulefiles/ccache/3.3.3
```

</div>

<div class='mac'>

```bash
brew install ccache
```

</div>

<div class='win'>

Windows support in ccache is not yet in the main ccache development line. Please install ccache from https://github.com/cristianadam/ccache/releases.

Just extract the archive and put the `ccache.exe` into the `PATH`.

</div>

You may want to change the cache directory (environment variable `CCACHE_DIR`) or increase the cache size (e.g. `ccache -M 10G` or `CCACHE_MAXSIZE`). See the [ccache docs](https://ccache.dev/manual/4.3.html#_configuration) for configuration instructions.

You can check cache hit statistics with `ccache -s`.

To disable caching:

```bash
cmake . -DOGS_DISABLE_COMPILER_CACHE=ON
```
