Path to the surface mesh vtu file.

The surface mesh must contain a nodal integer-valued field (unsigned 64bit)
named \c OriginalSubsurfaceNodeIDs. That field establishes the mapping between
the nodes of the surface mesh to some notes in the bulk mesh.
\warning It is not checked if the surface mesh and the bulk mesh correspond to each
other; in particular it is not checked if surface and bulk nodes coincide and if
surface elements coincide with the faces of bulk elements.
