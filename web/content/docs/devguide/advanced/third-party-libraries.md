+++
date = "2018-02-26T11:00:13+01:00"
title = "Third-party libraries"
author = "Lars Bilke"
weight = 1039

[menu]
  [menu.devguide]
    parent = "advanced"
+++

::: {.note}

### <i class="far fa-exclamation-triangle"></i> Attention

We strongly recommend to simply use [Conan]({{< ref "conan.pandoc" >}}) for handling required third-party libraries.
:::

## Introduction

The following libraries are required:

- [Boost](http://www.boost.org/) >= {{< dataFile "versions.minimum_version.boost" >}}
- [Eigen](http://eigen.tuxfamily.org/) >= {{< dataFile "versions.minimum_version.eigen" >}}
- [VTK](https://www.vtk.org/) >= {{< dataFile "versions.minimum_version.vtk" >}}

If you want to build the Data Explorer you need these too:

- [Qt](https://www.qt.io/) >= {{< dataFile "versions.minimum_version.qt" >}}
- [Shapelib](http://shapelib.maptools.org) >= 1.3.0
- [libgeotiff](https://trac.osgeo.org/geotiff/) >= 1.4.2

Please refer to the library documentation on how to build. If you struggle consider using [Conan]({{< ref "conan.pandoc" >}})! Once built CMake tries its best to find the libraries.
