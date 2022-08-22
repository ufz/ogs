+++
date = "2020-01-17T5:56:57+01:00"
title = "Extract Boundary"
author = "Thomas Fischer"
+++

## General

The tool extracts either lines in case of a 2D bulk mesh as input or
quads/triangles in case of a 3D bulk mesh as input. The input mesh can be given
either in the VTU or MSH format. Since the algorithm uses the element surface
normals a correct node ordering of the element nodes is required.

## Usage

```bash
ExtractBoundary -i [<file name of input mesh>] [-o <file name of output mesh>]
```

New data arrays containing the original IDs for nodes, elements, and faces from the source mesh are added to the extracted surface mesh. These data arrays are required for specifying source terms or setting boundary conditions during a simulation run of OpenGeoSys.

## Examples

### Extract the boundary from a quad mesh

`ExtractBoundary -i square_1x1_quad.vtu -o square_1x1_quad_border.vtu`

![The square mesh consists of 16 cells/elements.](ExtractBoundary_square_1x1_quad_border.png "The square mesh consists of 16 cells/elements. The numbers in the cells are the cell IDs. The generated boundary grid consists of the somewhat thicker and colored line elements.")

### Extract the boundary from a triangular mesh

`ExtractBoundary -i square_1x1_tri.vtu -o square_1x1_tri_border.vtu`

![The square mesh consists of 32 triangle shaped cells.](ExtractBoundary_square_1x1_tri_border.png "The square mesh consists of 32 triangle shaped cells. The numbers in the triangle are the cell IDs. The generated boundary grid consists of the somewhat thicker and colored line elements.")
