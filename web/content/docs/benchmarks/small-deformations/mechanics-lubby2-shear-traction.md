+++
project = "Mechanics/Burgers/cube_1e3.prj"
author = "Dmitri Naumov"
date = "2017-02-15T12:58:36+01:00"
title = "Lubby2; Creep example"
weight = 112

[menu]
  [menu.benchmarks]
    parent = "small-deformations"

+++

{{< project-link >}}

## Problem description

We solve a non-linear small deformation problem on a cube with shear traction on the top boundary face. The 3D problem is setup identical to the corresponding 2D problem.

See the TODO: {asset:795:link}-PDF for detailed problem description.

## Results and evaluation

Result showing the displacement field and distortion relative to the initial configuration:

{{< img src="../lubby2.png" >}}

Displacement of the top surface in the direction of the shear traction over time showing the elastic and viscous deformations (creep):

{{< img src="../lubby2_creep_over_time.png" >}}
