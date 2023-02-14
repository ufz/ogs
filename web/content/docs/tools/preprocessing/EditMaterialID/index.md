+++
date = ""
title = "EditMaterialID"
author = "Julian Heinze"
+++

## Description
editMaterialID is a tool to edit the material ID of mesh elements in three different ways. 
It is possible to:
- *replace* (-r, --replace) a current ID with a new one,
- *compress* (-c, --condense) the list of material IDs to their smallest possible values, e.g. [0,2,15,23,47] becomes [0,1,2,3,4],
- *specify* (-s, --specify) a new ID for all mesh elements of a certain group. 
The group refers to a certain element type like points, lines, quads, tetrahedra, hexahedra, triangles, prisms or pyramids.

## Usage
```bash
USAGE: 
   editMaterialID  {-r|-c|-s} [-e <point|line|quad|hex|tri|tet|pris|pyra>] 
                   [-n <number>][-m <number>] ... -o <file name> -i <file name [--] [--version] [-h]


Where: 
   -r,  --replace
     (OR required)  replace material IDs
         -- OR --
   -c,  --condense
     (OR required)  condense material IDs
         -- OR --
   -s,  --specify
     (OR required)  specify material IDs by element types (-e)


   -e <point|line|quad|hex|tri|tet|pris|pyra>,  --element-type <point|line
      |quad|hex|tri|tet|pris|pyra>
     element type

   -n <number>,  --new-material-id <number>
     new material id

   -m <number>,  --current-material-id <number>  (accepted multiple times)
     current material id to be replaced

   -o <file name>,  --mesh-output-file <file name>
     (required)  the name of the file the mesh will be written to

   -i <file name>,  --mesh-input-file <file name>
     (required)  the name of the file containing the input mesh

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.
```

## Example:
1. In this example we change all elements of material ID 1 to ID 14.

```bash
editMaterialID -i mesh_layerd.vtu -o mesh_newID.vtu -r -m 1 -n 14
```

**output:**
```bash
[2023-02-10 16:12:32.804] [ogs] [info] Mesh read: 2915 nodes, 4609 elements.
[2023-02-10 16:12:32.804] [ogs] [info] Replacing material ID...
[2023-02-10 16:12:32.804] [ogs] [info] 1 -> 14
```

2. In this example we compress/condense the list of material IDs. 
The input mesh has a list of material IDs of [0,4,5,6,14,42,53]. Then we call the --condense option (-c) to compress this list to [0,1,2,3,4,5,6]. 

```bash
editMaterialID -i mesh_newID.vtu -o mesh_newID_c.vtu -c
```

<p align='center'>
 <img src = newID.png width = "60%" height = "60%"> <img src = condense.png width = "60%" height = "60%"> 
</p>
<p align = "center">
Fig.1 The upper image shows a mesh with material IDs [0,4,5,6,14,42,53]. The lower one the same mesh with material IDs [0,1,2,3,4,5,6] after applying --condense.
 </p>

