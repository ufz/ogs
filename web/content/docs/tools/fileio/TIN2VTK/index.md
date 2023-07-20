+++
date = ""
title = "TIN2VTK"
author = "Julian Heinze"
+++

## Description

This tool converts datasets in TIN-format, usually used in geographic information systems (GIS), into VTK-format.
The TIN-format stores triangular irregular networks which can be considered a subclass of triangulated 2D meshes.
VTK stores such irregular networks as unstructured grids in `*.vtu`-files.

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

## Example

```bash
TIN2VTK -i input.tin -o output.vtu
```
