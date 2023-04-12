+++
date = ""
title = "AddFaultToVoxelGrid"
author = "Julian Heinze"
+++

## Description

This tool marks all elements in a voxel-grid (i.e. a structured hex grid, for instance created with [Layers2Grid](../../meshing/Layers2Grid/index.md) or [Vtu2Grid](../../meshing/vtu2grid/index.md)) that are intersected by a triangulated 2D mesh representing a fault or some other significant structure.
The material group for those intersected elements can be explicitly specified.
Otherwise the largest existing MaterialID will be increased by one and defined as MaterialID for those elements.

## Usage

```bash
AddFaultToVoxelGrid  -i <string> -f <string> -o <string>
                    [-m <non-negative integer>] [--] [--version] [-h]


Where:

   -i <string>,  --input <string>
     (required)  name of the input file list containing the paths the all
     input layers in correct order from top to bottom

   -f <string>,  --fault <string>
     (required)  name of mesh representing fault (*.vtu)

   -o <string>,  --output <string>
     (required)  name of output mesh (*.vtu)

   -m <non-negative integer>,  --material <non-negative integer>
     material id for cells intersected by fault

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

## Example

```bash
AddFaultToVoxelGrid -i grid.vtu -f fault.vtu -o grid_fault.vtu
```

<p align='center'>
 <img src = fault.png width = "100%" height = "100%">
</p>
<p align = "center">
Fig.1 The left figure shows the input grid with the intersecting triangulated 2D mesh, which is highlighted in red. The center and right figure show the output voxelgrid.
 </p>
