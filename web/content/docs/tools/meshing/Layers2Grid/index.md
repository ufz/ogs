+++
date = ""
title = "Layers2Grid"
author = "Julian Heinze"
+++

## Description
Reads a list of 2D unstructured mesh layers and samples them onto a structured grid of the same extent. 
The resulting mesh is referred to as voxelgrid.
Voxel sizes are defines by x/y/z-parameters. 
Note that a large cube size may result in an undersampling of the original structure.
For equilateral cubes, only the x-parameter needs to be set.


## Usage
```bash
Layers2Grid  -i <string> -o <string> -x <floating point number> 
            [-y <floating point number>] [-z <floating point number>]
            [-d] [--] [--version] [-h]


Where: 

   -i <string>,  --input <string>
     (required)  name of the input file list containing the paths the all
     input layers in correct order from top to bottom

   -o <string>,  --output <string>
     (required)  name of output mesh (*.vtu)

   -x <floating point number>,  --cellsize-x <floating point number>
     (required)  edge length of cubes in x-direction (longitude) or all
     directions, if y and z are not set

   -y <floating point number>,  --cellsize-y <floating point number>
     edge length of cubes in y-direction (latitude)

   -z <floating point number>,  --cellsize-z <floating point number>
     edge length of cubes in z-direction (depth)

   -d,  --dilate
     assign mat IDs based on single nodes instead of a majority of nodes,
     which can result in a slightly increased voxel grid extent

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

## Example:
In this example we will present how to create a voxel grid (output) from a list of layers in VTK-format (input). 
```bash
Layers2Grid -i layers.txt -o layers2grid.vtu -x 100 -y 200 -z 50
```
The list must contain the paths of the input layers as well as their names.
The voxelgrid is build in the same order as the names are listed from top to bottom. 

Example list of 10 layers:
```bash
path/to/layers/00_layer.vtu
path/to/layers/01_layer.vtu
path/to/layers/02_layer.vtu
path/to/layers/03_layer.vtu
path/to/layers/04_layer.vtu
path/to/layers/05_layer.vtu
path/to/layers/06_layer.vtu
path/to/layers/07_layer.vtu
path/to/layers/08_layer.vtu
path/to/layers/09_layer.vtu
path/to/layers/10_layer.vtu
```

<p align='center'>
 <img src = layers2grid.png width = "80%" height = "60%">
</p>
<p align = "center">
Fig.1 Voxelgrid composed of elements of size x=100, y=200 and z=50.
</p>

