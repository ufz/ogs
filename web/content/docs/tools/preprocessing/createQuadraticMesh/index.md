+++
date = ""
title = "createQuadraticMesh"
author = "Julian Heinze"
+++

## Description
This tool converts a linear mesh to a mesh of quadratic order.
The order of a mesh is equal to the highest order of elements within the mesh.
The quadratic order mesh allows shape functions of quadratic order to be considered in the finite element methods of the simulation.
More detailed information about the differences between a linear and quadratic mesh elements can be found in the [VTK book by kitware](https://kitware.github.io/vtk-examples/site/VTKBook/05Chapter5/).

## Usage
```bash
USAGE: 
   createQuadraticMesh  [-c] -o <string> -i <string> [--] [--version]
                            [-h]


Where: 
   -c,  --add-centre-node
     add centre node

   -o <string>,  --output-mesh-file <string>
     (required)  output mesh file

   -i <string>,  --input-mesh-file <string>
     (required)  input mesh file

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

## Example:

```bash
createQuadraticMesh -i input_mesh.vtu -o quadratic_mesh.vtu
```