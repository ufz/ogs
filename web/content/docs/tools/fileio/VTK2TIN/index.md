+++
date = ""
title = "VTK2TIN"
author = "Julian Heinze"
+++

## Description
This tool converts VTK unstructured grids (*.vtu) into TIN-format usable in geographic information systems (GIS).
The TIN-format stores triangular irregular networks which can be considered a subclass of triangulated 2D meshes. 
The vtu-format can store a large variety of unstructured mesh types, but only 2D triangle meshes can be converted using this tool.
## Usage
```bash
USAGE: 
   TIN2VTK  -o <string> -i <string> [--] [--version] [-h]


Where: 
   -o <string>,  --output-vtu-file <string>
     (required)  the name of the file the mesh will be written to

   -i <string>,  --input-tin-file <string>
     (required)  the name of the file containing the input TIN

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

## Example:

```bash
VTK2TIN -i input.vtu -o output.tin
```