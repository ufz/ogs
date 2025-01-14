+++
date = "2018-02-27T11:00:13+01:00"
title = "Complete workflow for simulating a geological model with time and space dependent boundary conditions (advancing glacier)"
author = "Dominik Kern"
weight = 1
image = "advancing_glacier_thumbnail.png"
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

## Supplementary material

<!-- vale off -->

* [mesh_basin.msh](https://ogsstorage.blob.core.windows.net/web/tutorials/advancing-glacier/mesh_basin.msh)
* [mesh_basin.py](mesh_basin.py)
* [OGSinput_basin.prj](OGSinput_basin.prj)
* [timeBCs_glacier.py](timeBCs_glacier.py)
* [glacierclass.py](glacierclass.py)
* [history.sh](history.sh)
* [msh2vtu](https://gitlab.opengeosys.org/ogs/tools/ogstools) *-now integrated in ogstools*
* *Transient hydrodynamics within intercratonic sedimentary basins during glacial cycles*, V. F. Bense  M. A. Person, DOI: [10.1029/2007JF000969](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007JF000969)

***NOTE:***  Due to ongoing development, msh2vtu.py was integrated in the ogs-repository for [ogstools](https://gitlab.opengeosys.org/ogs/tools/ogstools).
In order to work with the tutorial and apply msh2vtu on mesh_basin.msh, ogstools need to be installed at first.

ogstools can be installed from PyPI using pip.
It is recommended to set up a virtual environment beforehand:

```bash
python -m venv .venv
source .venv/bin/activate
pip install ogstools
```

More detailed information about the usages and development of ogstools can be found in [ogstools Documentation](https://ogstools.opengeosys.org/stable/index.html).
After installation the tool can be applied as described in the first video-tutorial.
To apply msh2vtu on the mesh_basin.msh the command should look like:

```bash
msh2vtu --ogs --rdcd mesh_basin.msh
```
