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

See [this PDF](../LIE_small_deformation.pdf) for detailed problem description.

The one-sided incompressibility constraint for fracture models is described in
[this PDF](../LIE_fracture_incompressibility.pdf).

## Results and evaluation

Result showing sliding of the upper part of the domain along the fracture:

{{< img src="../LIE_SD_m_result_uy.png" >}}


Same benchmark setup with plane strain conditions in 3D:
{{< img src="../single_joint_3D.png" >}}

Comparison with 2D setup yields identical results (up to numerical differences
in order of 1e-15); Resulting displacement on the left axis, and error on the
right:
{{< img src="../single_joint_3D_2D_results.png" >}}
