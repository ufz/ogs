+++
author = "[TODO]"
project = "[TODO].prl"

[menu]
  [menu.benchmarks]
    parent = "[TODO]"
    weight = 1
+++

{{< project-link >}}

## Problem description

## Input files

The main project file is `[TODO]`. It describes the processes to be solved and the related process variables together with their initial and boundary conditions. It also references the mesh and geometrical objects defined on the mesh.

As of now a small portion of possible inputs is implemented; one can change:
 - the mesh file
 - the geometry file
 - introduce more/different Dirichlet boundary conditions (different geometry or values)

The geometries used to specify the boundary conditions are given in the `[TODO]` file.

The input mesh `[TODO]` is stored in the VTK file format and can be directly visualized in Paraview for example.

## Running simulation

To start the simulation (after successful compilation) run:
```bash
$ ogs [TODO]
```

## Results and evaluation
