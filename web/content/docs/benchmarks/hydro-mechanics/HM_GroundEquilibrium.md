+++
author = "Christian Silbermann, Thomas Nagel"
project = "HydroMechanics/GroundEquilibrium"
date = "2021-12-14T14:39:39+01:00"
title = "Ground equilibrium"

[menu]

  [menu.benchmarks]
    parent = "hydro-mechanics"

+++

## Test cases

Four versions of the test are presented:
 - {{< data-link "Gravity load" "HydroMechanics/GroundEquilibrium/simHM_ground.prj" >}},
 - {{< data-link "Gravity load with Python BCs" "HydroMechanics/GroundEquilibrium/simHM_ground_python.prj" >}},
 - {{< data-link "Gravity load and indentation" "HydroMechanics/GroundEquilibrium/simHM_ground_quadBCu.prj" >}}, and
 - {{< data-link "Gravity load and indentation with Python BCs" "HydroMechanics/GroundEquilibrium/simHM_ground_quadBCu_python.prj" >}}.

## Problem description

We assume a two-dimensional piece of ground (linear elastic soil saturated with water) and consider two cases: The static equilibrium state under gravity load, and the temporal evolution under some additional indentation from the top. The corresponding boundary conditions are either coded directly in the prj-file or provided by *pythonBCsOGS.py* using the Python interface.
See [this PDF](../HM_GroundEquilibrium.pdf) for the detailed description.
