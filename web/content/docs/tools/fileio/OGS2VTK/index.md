+++
date = ""
title = "OGS2VTK"
author = "Julian Heinze"
+++

## Description

OGS2VTK is a tool to convert OGS-mesh-files to VTK-files.
It can be applied to format OGS-5 legacy mesh files or visualize them in ParaView.

## Usage

```bash
USAGE:
   OGS2VTK  [--ascii_output] -o <file name of output mesh> -i <file name of input mesh> [--] [--version] [-h]

Where:

   --ascii_output
     Write VTU output in ASCII format. Due to possible rounding the ascii
     output could result in lower accuracy.

   -o <file name of output mesh>,  --mesh-output-file <file name of output mesh>
     (required)  the name of the file the mesh will be written to

   -i <file name of input mesh>,  --mesh-input-file <file name of input mesh>
     (required)  the name of the file containing the input mesh

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

## Example

Converting a `.msh`-file to a `.vtu`-file.

```bash
OGS2VTK -i example.msh -o example.vtu
```
