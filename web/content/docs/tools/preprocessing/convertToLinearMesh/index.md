+++
date = ""
title = "convertToLinearMesh"
author = "Julian Heinze"
+++

## Description
This tool allows to convert a non-linear mesh to a linear mesh.
All additional nodes from higher order elements in the non-linear mesh are removed when writing the linear mesh, so the meshes can be simplified.
This way, simulation results from non-linear approaches to linear ones can be compared. 

## Usage
```bash
USAGE: 
   convertToLinearMesh  -o <string> -i <string> [--] [--version] [-h]


Where: 
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
convertToLinearMesh -i non-linear.vtu -o linear.vtu
```