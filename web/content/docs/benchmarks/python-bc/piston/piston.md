+++
project = "Mechanics/Linear/PythonPiston/piston.prj"
author = "Christoph Lehmann"
date = "2018-08-16T12:18:00+02:00"
title = "Ideal Gas Compressed by an Elastic Piston"
weight = 2

[menu]
  [menu.benchmarks]
    parent = "python-bc"

+++

{{< data-link >}}

## Motivation of this test case

The aim of this test is:

* to show that it is possible to prescribe BCs only on parts of a given geometry
* to assert that the computation of fluxes that depend nonlinearly on the
  primary variables is implemented correctly.

## Details

A chamber filled with ideal gas is sealed tightly with a movable, elastic
piston. The position of the piston is varied between different load steps.
Friction between the piston and the chamber wall is neglected.
For simplicitly, also initially the elastic piston is in an unstressed state.

{{<img src="../sketch-piston.png" >}}


## Results

{{<img src="../load-steps.png" >}}

The figure above shows that the piston is being compressed
($y$ displacement has larger negative values at the top)
by the forces acting on it.
The initial position of the top part of the piston is indicated as a wireframe.

{{<img src="../pressure-displacement.png" >}}

The plot shows that the relation between the stress in the piston and its
displacement coincides with the pressure-volume relation of the chamber.
That indicates that the Jacobian computation inside the Python BC classes was
correctly implemented.
