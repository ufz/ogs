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

To disable caching with ccache:

```bash
cmake . -DOGS_DISABLE_CCACHE=ON
```

## Usage on eve

Just load the module:

```bash
module load /global/apps/modulefiles/ccache/3.3.3
```
