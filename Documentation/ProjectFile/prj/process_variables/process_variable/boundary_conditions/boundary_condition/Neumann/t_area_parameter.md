Typically boundary conditions are set on domains with dimension exactly one
lower than the bulk mesh dimension. In some cases the boundary condition should
be set on entities that have lower dimension than one lower. For instance
boundary conditions on points in 2d or 3d domains or boundary conditions on
lines in 3d domains.

In order not to integrate over a null set the area tag can be used to specify
the area the boundary condition should be effective.
