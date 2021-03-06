+++
date = "2021-03-06T00:00:00+01:00"
title = "GMSH2OGS"
author = "Dmitri Naumov"

[menu]
  [menu.tools]
    parent = "Data Import/Export"
+++

## Introduction

Meshes generated with Gmsh can be converted to OGS meshes, in particular to VTU
file format.
Currently, the only supported Gmsh format is 2.1, which can be specified for
example with the command line flag `-format msh2` when executing `gmsh`.

## Usage

Gmsh writes all elements (unless specified otherwise) including the lower
dimensional elements, usually model boundaries.
Whereas for OGS we separate the so called "bulk" mesh and the boundary (or
subdomain) meshes.
`GMSH2OGS` options `-e` (`--exclude-lines`) and `-b` (`--boundaries`) control
whether the line elements are removed or written in separate files.

## Example

A 2D-mesh with two physical groups (MaterialIDs) was generated with Gmsh.
It also includes the outer boundary and inner boundary line elements.

### Input files and results

All files are stored in
[Tests/Data/MeshLib/](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/MeshLib):
 - the Gmsh generated mesh
[A2-gmsh.msh](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/MeshLib/A2-gmsh.msh),
 - and the corresponding result files
[A2.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/MeshLib/A2.vtu),
and the boundary meshes
[A2_0.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/MeshLib/A2_0.vtu) to
[A2_7.vtu](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/MeshLib/A2_7.vtu).

### Running GMSH2OGS

`GMSH2OGS` creates different meshes from the same input depending on the `-e`
and `-b` flags; the results are shown (z-translated) in the figure below.

```bash
GMSH2OGS -i A2-gmsh.msh -o A2.vtu
```

 - Conversion without any additional flags yields a single VTU file with
   outer and inner boundaries (line elements) included, shown in the very
   bottom.
 - Adding only a `-e` flag results in a single VTU file, now without the line
   elements.
 - Adding only a `-b` flag gives a single VTU file as in the first case *and*
   additional eight files, each containing a line-element mesh corresponding to
   the physical groups;
   these are the white lines in the figure, shown again z-translated.
 - Finally specifying both flags (`-e` and `-b`) produces a single VTU file
   without the line elements and additional eight boundary (subdomain) files.

![GMSG2OGS-e-b](./extract_boundary.png#two-third "GMSH2OGS meshes for -e and -b
command line flags.")

