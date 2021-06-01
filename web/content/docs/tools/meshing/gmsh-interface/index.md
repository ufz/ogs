+++
date = "2021-05-31T11:00:13+01:00"
title = "geometryToGmshGeo"
author = "Thomas Fischer"
aliases = [ "/docs/tools/meshing/gmsh-interface" ]

[menu]
  [menu.tools]
    parent = "meshing"
+++

## Introduction

The tool `geometryToGmshGeo` takes OGS geometries (gml files) and creates a Gmsh
geo file. The user can specify several command line arguments to influence the
creation procedure. A list of the arguments as well as a description of each
argument can be obtained with the `--help` option.

## Usage

```bash
geometryToGmshGeo --input <input file name> --output <output file>
```
The `--input` argument is accepted multiple times.

## Simple examples

### First Example: Simple Geometry

First, the Gmsh geometry file is create from the gml file [square_1x1.gml](square_1x1.gml):
```bash
geometryToGmshGeo -i square_1x1.gml -o /tmp/square_1x1.geo --mesh_density_scaling_at_points 0.05
```

Then, the Gmsh geometry can be meshed:
```bash
gmsh /tmp/square_1x1.geo -2 -algo meshadapt -format msh22 -o /tmp/square_1x1.msh
```

![geometry](square_1x1.gml.png#OneThird "1x1 square geometry")
![coarse mesh](square_1x1_adaptive_point_density_0.5.png#OneThird "simple coarse mesh (density scaling 0.5)")
![fine mesh](square_1x1_adaptive_point_density_0.05.png#OneThird "simple fine mesh (density scaling 0.05)")
![even finer mesh](square_1x1_adaptive_point_density_0.005.png "even finer mesh (density scaling 0.005)")

### Second Example: Simple Geometry with Additional Geometrical Information

First, the Gmsh geometry file is create from the gml files
 - [square_1x1.gml](square_1x1.gml)
 - [square_0.15_0.25x0.15_0.25.gml](square_0.15_0.25x0.15_0.25.gml)
 - [square_0.45_0.55x0.45_0.55.gml](square_0.45_0.55x0.45_0.55.gml)

 ```bash
geometryToGmshGeo -i square_1x1.gml
                   -i square_0.15_0.25x0.15_0.25.gml
                   -i square_0.45_0.55x0.45_0.55.gml
                   -o /tmp/square_1x1.geo --mesh_density_scaling_at_points 0.5
```

Then, the Gmsh geometry can be meshed:
```bash
gmsh /tmp/square_1x1.geo -2 -algo meshadapt -format msh22 -o /tmp/square_1x1.msh
```

![geometry](square_1x1_plus_subgeometries.png#OneThird "1x1 square geometry and sub geometries")
![coarse mesh](square_1x1_plus_subgeometries_adaptive_point_density_0.5.png#OneThird "coarse mesh (density scaling 0.5)")
![fine mesh](square_1x1_plus_subgeometries_adaptive_point_density_0.05.png#OneThird "fine mesh (density scaling 0.05)")
