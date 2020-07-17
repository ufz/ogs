+++
date = "2018-10-02T15:56:57+01:00"
title = "Identify subdomains in bulk mesh."
author = "Dmitri Naumov"

[menu]
  [menu.tools]
    parent = "meshing-submeshes"
+++

## General

`identifySubdomains` is a creation and check tool for the `bulk_node_ids` and
`bulk_element_ids` mesh properties. These properties are needed for the boundary
conditions and source terms application, defined on the subdomains of a "bulk"
mesh.

## Example

Given a "bulk" mesh (Tests/Data/Mechanics/Linear/disc_with_hole.vtu) and a
[quater cirle mesh](quater_circle.vtu) extracted manually we want to use the
quater circle mesh for heterogeneous boundary condition. OGS requires two
mappings into the "bulk" mesh, one for the nodes and one for the elements.

![The figure shows a part of the "bulk" mesh with boundary element numbers, and
the quater circle mesh shown as white line with green
points.](disc_with_hole_and_bondary.png){width=50% .m-auto}

To create this mappings we run

```bash
identifySubdomains -m Tests/Data/Mechanics/Linear/disc_with_hole.vtu -s 1e-6 -o
new_ -- quater_circle.vtu
```

The tool will first try to find all unique nodes in the "bulk" mesh using search
radius 1e-6, and create the `bulk_node_ids` mapping upon success. Then the
`bulk_element_ids` mapping is created by finding a unique element containing all
the nodes of the subdomain element. The output file
[`new_quater_circle.vtu`](new_quater_cirle.vtu) will now contain both
mappings and is prepared for usage as a boundary condition mesh.

## Notes

- The double dash is separating the subdomain meshes, so the input can have
   multiple subdomain inputs:

   ```bash
   identifySubdomains -m bulk.vtu -- bc1.vtu bc2.vtu source_term_meshes*.vtu
   ```

- The output prefix `-o` can contain a path too and is relative to the current
   working directory.
