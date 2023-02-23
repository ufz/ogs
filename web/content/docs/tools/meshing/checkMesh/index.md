+++
date = ""
title = "checkMesh"
author = "Julian Heinze"
weight = 1
+++

## Description
This tool can be used to check the validity of the input mesh and optionally output mesh properties (specifically MaterialIDs).
## Usage
```bash
USAGE: 
   checkMesh  [-p] [-v] [--] [--version] [-h] <string>


Where: 
   -p,  --print_properties
     print properties stored in the mesh

   -v,  --validation
     validate the mesh

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <string>
     (required)  input mesh file

```

## Example:
In this example we use a 3D mesh from a previous examples ([createLayeredMeshfromRaster](/docs/tools/preprocessing/createLayeredMeshFromRasters/index.md)).

```bash
checkMesh mesh_layerd.vtu -p -v
```

 **general output:**
 At first there are some general information given.
 ```bash
[2023-02-10 11:44:56.551] [ogs] [info] Memory size: 1 MiB
[2023-02-10 11:44:56.551] [ogs] [info] Time for reading: 0.056879 s
[2023-02-10 11:44:56.552] [ogs] [info] Axis aligned bounding box: 	x [394412, 395388) (extent 976.525)
	y [5.82349e+06, 5.82454e+06) (extent 1052.17)
	z [-120, 57.5427) (extent 177.543)
[2023-02-10 11:44:56.552] [ogs] [info] Min/max edge lengths: [0.0123914, 97.0344]
[2023-02-10 11:44:56.552] [ogs] [info] Number of elements in the mesh:
[2023-02-10 11:44:56.552] [ogs] [info] 	Tetrahedrons: 88
[2023-02-10 11:44:56.552] [ogs] [info] 	Pyramids: 95
[2023-02-10 11:44:56.552] [ogs] [info] 	Prisms: 4426
 ```
 **-p output:**
 The properties output refers to material properties. The bounds of the materialID vector are [0,6].
 ```bash
[2023-02-10 11:39:56.092] [ogs] [info] There are 1 properties in the mesh:
[2023-02-10 11:39:56.092] [ogs] [info] 	MaterialIDs: (4609 values) [0, 6]
 ```
 **-v output:**
 The quality of the underlying mesh is checked and no errors are found.
  ```bash
[2023-02-10 11:57:51.477] [ogs] [info] Mesh Quality Control:
[2023-02-10 11:57:51.477] [ogs] [info] Looking for unused nodes...
[2023-02-10 11:57:51.478] [ogs] [info] Found 0 potentially collapsible nodes.
[2023-02-10 11:57:51.479] [ogs] [info] Testing mesh element geometry:
[2023-02-10 11:57:51.483] [ogs] [info] No errors found.
[2023-02-10 11:57:51.483] [ogs] [info] No elements found with zero volume.

[2023-02-10 11:57:51.483] [ogs] [info] No elements found with non coplanar nodes.

[2023-02-10 11:57:51.483] [ogs] [info] No elements found with non-convex geometry.

[2023-02-10 11:57:51.483] [ogs] [info] No elements found with wrong node order.

[2023-02-10 11:57:51.486] [ogs] [info] 1 property vectors copied, 0 vectors skipped.
[2023-02-10 11:57:51.486] [ogs] [info] No holes found within the mesh.
 ```

[//]: # (Note: Inconsistent usage of input mesh. Most other tools take -i as input.)