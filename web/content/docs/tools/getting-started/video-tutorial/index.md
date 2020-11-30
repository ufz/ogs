+++
date = "2018-02-27T11:00:13+01:00"
title = "Video Tutorial"
author = "Dominik Kern"
weight = 105

[menu.docs]
name = "Tools & Workflows"
identifier = "tools"
weight = 4
post = "Helpful tools for pre- and postprocessing as well as complete model setup workflows."

[menu.tools]
parent = "getting-started"
+++

In this [three-part tutorial](https://www.youtube.com/watch?v=BULunRJQRJ0&list=PLU_clTnZqNAeOXENl79kQwn0pgHGittX1) we look at the mechanical deformation of a sedimentary basin under an advancing glacier. In detail we show you a possibility how to

* create a mesh with [Gmsh](http://gmsh.info/),
* implement space- and time-dependent boundary conditions,
* write an OGS input file,
* visualize the results with [ParaView](https://www.paraview.org/).

As a common programming language we use [Python](https://www.python.org).

## Part 1: Preprocessing

{{< youtube BULunRJQRJ0 >}}

## Part 2: Solving

{{< youtube GL5sugIyHEk >}}

## Part 3: Postprocessing

{{< youtube bkmubABAA_s >}}


## Supplementary material:

* [mesh_basin.msh](mesh_basin.msh)
* [mesh_basin.py](mesh_basin.py)
* [OGSinput_basin.prj](OGSinput_basin.prj)
* [timeBCs_glacier.py](timeBCs_glacier.py)
* [glacierclass.py](glacierclass.py)
* [history.sh](history.sh)
* [msh2vtu-project on GitHub](https://github.com/dominik-kern/msh2vtu)
* *Transient hydrodynamics within intercratonic sedimentary basins during glacial cycles*, V. F. Bense  M. A. Person, DOI: [10.1029/2007JF000969](https://doi.org/10.1029/2007JF000969)
