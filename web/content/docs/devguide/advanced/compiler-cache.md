+++
date = "2018-02-26T11:00:13+01:00"
title = "Using a compiler cache"
author = "Lars Bilke"
weight = 1038

aliases = ["/docs/devguide/advanced/using-ccache"]

[menu]
  [menu.devguide]
    parent = "advanced"
+++

## Introduction

A compiler cache speeds up compilation times by caching object files and reusing them on subsequent builds. This even works for complete rebuilds (i.e. deleting the full build-directory). A compiler cache is automatically used when it is found by CMake.

<div class='linux'>

For Linux and macOS you can use [ccache](https://ccache.samba.org). You may want to change the cache directory (environment variable `CCACHE_DIR`) or increase the cache size (`CCACHE_MAXSIZE`). See the [ccache docs](https://ccache.dev/manual/4.2.html#_configuration) for configuration instructions.

You can check cache hits with `ccache -s`.

## Usage on eve

Just load the module:

```bash
module load /global/apps/modulefiles/ccache/3.3.3
```

</div>

<div class='mac'>

See Linux-tab!

</div>

<div class='win'>

On Windows you can use [clcache](https://github.com/frerich/clcache). Install it with [Chocolatey](https://github.com/frerich/clcache/wiki/Installation#installing-an-executable-file-via-chocolatey):

```bash
curl -L https://github.com/frerich/clcache/releases/download/v4.2.0/clcache.4.2.0.nupkg --output clcache.4.2.0.nupkg
choco install clcache --source=.
```

You may want to increase the cache size (to 10 GB in this case):

```bash
clcache -M 10000000000
```

You can check cache hits with `clcache -s`.

</div>

To disable caching:

```bash
cmake . -DOGS_DISABLE_COMPILER_CACHE=ON
```
