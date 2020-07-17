+++
date = "2020-01-17T5:56:57+01:00"
title = "Extract Boundary"
author = "Thomas Fischer"

[menu]
  [menu.tools]
    parent = "meshing-submeshes"
+++

## General

The tool extracts either lines in case of a 2D bulk mesh as input or
quads/triangles in case of a 3D bulk mesh as input. The input mesh can be given
either in the vtu or msh format. Since the algorithm uses the element surface
normals a correct node ordering of the element nodes is required.

## Usage

```bash
ExtractBoundary -i [<file name of input mesh>] [-o <file name of output mesh>]
    [--face-property-name <string>]
    [--element-property-name <string>]
    [--node-property-name <string>]
```

The data arrays added to the boundary mesh by using the options
`--face-property-name` (default value 'bulk_face_ids'),
`--element-property-name` (default value 'bulk_element_ids'),
and `--node-property-name` (default value 'bulk_node_ids')
are used in other tools (for instance in
[ComputeNodeAreasFromSurfaceMesh]({{< ref "compute-node-areas-from-surface-mesh" >}}))
and are required for flux calculations during a simulation run of OpenGeoSys.

## Examples

### Extract the boundary from a quad mesh

`ExtractBoundary -i square_1x1_quad.vtu -o square_1x1_quad_border.vtu`

![The square mesh consists of 16 cells/elements. The numbers
in the cells are the cell IDs. The generated boundary grid consists of the
somewhat thicker and colored line elements.](ExtractBoundary_square_1x1_quad_border.png){.m-auto}

### Extract the boundary from a tri mesh

`ExtractBoundary -i square_1x1_tri.vtu -o square_1x1_tri_border.vtu`

![The square mesh consists of 32 triangle shaped cells. The
numbers in the tri are the cell IDs. The generated boundary grid consists of the
somewhat thicker and colored line elements.
](ExtractBoundary_square_1x1_tri_border.png){.m-auto}
