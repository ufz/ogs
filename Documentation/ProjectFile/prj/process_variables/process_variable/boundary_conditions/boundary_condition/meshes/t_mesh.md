Name of the subdomain mesh where the boundary condition will be defined.

The subdomain mesh must contain a nodal integer-valued field (unsigned 64bit)
named `bulk_node_ids` and a cell field named `bulk_element_ids`.
These fields establish the mapping between the nodes of the surface mesh to the
nodes in the bulk mesh.

\warning It is not checked if the surface mesh and the bulk mesh
correspond to each other; in particular it is not checked if surface and bulk
nodes coincide and if surface elements coincide with the faces of bulk elements.
