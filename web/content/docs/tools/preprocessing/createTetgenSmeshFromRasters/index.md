+++
date = "2023-06-22"
title = "createTetgenSmeshFromRasters"
author = "Julian Heinze"
+++

## Description

createTetgenSmeshFromRasters is a tool for creating a TetGen- `*.smesh` file from a 2D input VTK mesh and one or more raster files.
The 2D input file is intended to define the upper surface of the resulting `*.smesh`, while the raster files are intended to describe the geologic layers below this surface.
The resulting `*.smesh` file is a *piecewise linear complex* (PLC) that describes a boundary representation for a layered 3D mesh.
This PLC serves as input to TetGen's *Tetrahedral Mesh Generator* to create a 3D mesh.
A more detailed explanation of the use and functionality of TetGen can be found in the manual.
Supported raster formats are ArcGIS ASCII raster (`*.asc`), Surfer grids (`*.grd`), and XYZ raster files (`*.xyz`).
Still, the list of raster files can be given in any text file format.

### Known limitations

Currently, only inputs consisting of line and triangle elements are supported, as mapping quads may result in invalid grid elements.

## Usage

```bash
createTetgenSmeshFromRasters  -i <file name> -o <file name> -r <file
                                     name> [-t <floating point number>]
                                     [--ascii_output] [--] [--version]
                                     [-h]

Where:

   -i <file name>,  --input-mesh-file <file name>
     (required)  The file name of the 2D input mesh.

   -o <file name>,  --output-mesh-file <file name>
     (required)  The file name of the resulting 3D mesh.

   -r <file name>,  --raster-list <file name>
     (required)  An ascii-file containing a list of raster files, starting
     from top (DEM) to bottom.

   -t <floating point number>,  --thickness <floating point number>
     The minimum thickness of a layer to be integrated at any given
     location.

   --ascii_output
     Write VTU output in ASCII format.

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

## Example

```bash
createLayeredMeshFromRasters -i mesh.vtu -o mesh_layered.smesh -r list_rasters.txt
```

The text file that contains the list of raster files, in this example it is called "list_rasters.txt", is specified as follows:

```bash
path/to/raster-file/layer_0.asc
path/to/raster-file/layer_1.asc
path/to/raster-file/layer_2.asc
path/to/raster-file/layer_3.asc
path/to/raster-file/layer_4.asc
path/to/raster-file/layer_5.asc
```

**Note:** If it is wanted to create a 3D mesh from the resulting `mesh_layered.smesh` file, it is recommended to use [TetGen](https://wias-berlin.de/software/tetgen/). The `tetgen`-binary can be installed with your package manager or with [Conda]:

```bash
conda config --add channels conda-forge
conda install tetgen
```

Then it could be called with e.g.:

```bash
tetgen -pkA mesh_layered.smesh
```

The **p** takes care of the tetrahedralization, **k** takes care of the generation of a VTK output, and **A** writes cell properties for regions to the output mesh.
