+++
weight = 114
project = "LIE/Mechanics/single_joint.prj"
author = "Norihiro Watanabe"
date = "2017-02-15T14:43:32+01:00"
title = "Linear; Single fracture"

[menu]
  [menu.benchmarks]
    parent = "small-deformations"

+++

{{< project-link >}}

## Problem description

We solve a linear elastic small deformation problem with a pre-existing fracture using the lower-dimensional interface element (LIE) approach.

See the [LIE_small_deformation.pdf](https://docs.opengeosys.org/assets/files/Documentation/Selected-Benchmarks/LIE_small_deformation.pdf) for detailed problem description.

## Results and evaluation

Result showing sliding of the upper part of the domain along the fracture:

{{< img src="../LIE_SD_m_result_uy.png" >}}
