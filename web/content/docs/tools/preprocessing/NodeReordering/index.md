+++
date = "2023-02-23"
title = "NodeReordering"
author = "Julian Heinze"
+++

## Description

This tool is to reorder nodes of a given mesh.
Mostly it is applied to reorder the nodes of a mesh from an older OGS version to have a mesh compatible with the current OGS version.
There are four different methods available, which define the reordering process:

- Method 0: Same as Method 1, but additionally verifying the Jacobian
determinant of the reordered elements.
- Method 1: Reversing order of nodes unless it is considered correct by
OGS-6, i.e. such that det(J) > 0, where J is the Jacobian of the
global-to-local coordinate transformation.
This is the default selection.
- Method 2: Fixing node ordering issues between VTK and OGS-6 (only applies
to prism elements).
- Method 3: Re-ordering of the mesh node vector such that all base nodes are
sorted before all nonlinear nodes.

With Methods 0 and 1 the volumes of the original and the reordered elements
are compared, and the tool aborts with a diagnostic if both the absolute and
the relative volume difference exceed 100 machine epsilons;
`--no_volume_check` disables that comparison.

## Usage

```bash
USAGE:

   NodeReordering  [-h] [--no_volume_check] [--version] [-l <none|error|
                   warn|info|debug|all>] [-m <0|1|2|3>] -i <INPUT_FILE> -o
                   <OUTPUT_FILE>


Where:

   --no_volume_check
     By default the volumes of original and reordered elements are compared
     if they are numerically equal, i.e., if the absolute difference or the
     relative volume difference is at most 100 machine epsilons. This
     switch disables the volume comparison.

   -m <0|1|2|3>,  --method <0|1|2|3>
     reordering method selection

   -o <OUTPUT_FILE>,  --output_mesh <OUTPUT_FILE>
     (required) Output (.vtu). The name of the output mesh file

   -i <INPUT_FILE>,  --input_mesh <INPUT_FILE>
     (required) Input (.vtu | .vtk | .msh). The name of the input mesh file

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

The following command is used to reorder a mesh according to Method 2, i.e.
to fix node ordering issues between VTK and OGS-6.

```bash
NodeReordering -i input.vtu -o output.vtu -m 2
```
