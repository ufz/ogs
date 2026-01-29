+++
date = "2023-02-23"
title = "NodeReordering"
author = "Julian Heinze"
+++

## Description

This tool is to reorder nodes of a given mesh.
Mostly it is applied to reorder the nodes of a mesh from an older OGS version to have a mesh compatible with the current OGS version.
There are four different methods available, which define the reordering process:

- Method 0: Reversing order of nodes for all elements.
- Method 1: Reversing order of nodes unless it's perceived correct by OGS-6 standards.
This is the default selection.
- Method 2: Fixing node ordering issues between VTK and OGS-6 (only applies to prism-elements).
- Method 3: Re-ordering of mesh node vector such that all base nodes are sorted before all nonlinear.nodes.

## Usage

```bash
USAGE:
   NodeReordering  -i <filename> -o <filename> [-m <0|1|2|3>] [--]
                       [--version] [-h]

Where:
   -i <filename>,  --input_mesh <filename>
     (required)  the name of the input mesh file

   -o <filename>,  --output_mesh <filename>
     (required)  the name of the output mesh file

   -m <0|1|2|3>,  --method <0|1|2|3>
     reordering method selection

   --no_volume_check
     By default the volumes of original and reordered elements are
     compared if they are numerically equal, i.e., relative volume
     difference is smaller than a threshold. This switch disables the
     volume comparison.

   -l <none|error|warn|info|debug|all>,  --log-level <none|error|warn|info|
      debug|all>
     the verbosity of logging messages

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

## Example

The following command is used to reorder a mesh according to the second method.

```bash
NodeReordering -i input.vtu -o output.vtu -m 2
```
