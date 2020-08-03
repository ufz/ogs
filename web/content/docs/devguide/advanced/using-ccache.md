+++
date = "2018-02-26T11:00:13+01:00"
title = "Using ccache"
author = "Lars Bilke"
weight = 1038

[menu]
  [menu.devguide]
    parent = "advanced"
+++

<div class='note'>

### <i class="far fa-exclamation-triangle"></i> GCC-like compilers only

Tested on GCC and Clang.
</div>

## Introduction

[ccache](https://ccache.samba.org) speeds up compilation times by caching object files and reusing them on subsequent builds. This even works for complete rebuilds (i.e. deleting the full build-directory). If ccache is found on the system it is automatically used.

## Configuration

Add the following line to your [ccache config file](https://ccache.samba.org/manual.html#_configuration) which is required for pre-compiled headers:

```bash
sloppiness = pch_defines,time_macros
```

See the [ccache docs](https://ccache.samba.org/manual.html#_configuration_settings) for other available options.

### ccache and Clang

Set the option `run_second_cpp = true` or `export CCACHE_CPP2=true` to suppress lots of [false positive compiler warnings](http://peter.eisentraut.org/blog/2014/12/01/ccache-and-clang-part-3/).

## Usage on eve

Just load the module:

```bash
module load /global/apps/modulefiles/ccache/3.3.3
```
