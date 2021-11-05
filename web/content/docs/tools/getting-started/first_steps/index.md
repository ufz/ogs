+++
date = "2019-03-14T11:00:13+01:00"
title = "First Steps"
author = "Joerg Buchwald"
weight = 101

[menu]
  [menu.tools]
    parent = "getting-started"
+++

## First Steps

This section describes a possible general pre-processing workflow using some of the tools provided by OGS.

To set up a model domain, the project file `*.prj` requires additional files for the mesh and geometry to apply the boundary conditions accordingly.

The bulk mesh must be provided in VTK's `*.vtu` format, whereas the boundary and source term domains can be provided either in OGS' internal geometry file format (filename extension: `*.gml`, not to be confused with the Geography Markup Language) or as a `*.vtu` file as well, containing a subdomain of the bulk mesh. Multiple `*.vtu` files can be provided in the project file using the `<meshes>` tag. We recommend the usage of the first method for simple 2D meshes with constant boundary conditions, whereas more complicated geometries and conditions might require the latter method. One general advantage in the utilization of `*.vtu` files is that they allow a definition of additional field variables at each mesh node/element in order to implement spatially varying boundary conditions in a similar manner as, e.g., defining inhomogeneous material properties in the bulk mesh.

The finite-element mesh must be created using external mesh generators. Simple meshes can be created using the [generateStructuredMesh]({{< ref "structured-mesh-generation" >}}) tool. More complicated geometries can be generated and meshed using, e.g., [SALOME Platform](https://www.salome-platform.org) or [GMSH](http://gmsh.info).

To use a mesh created by SALOME it must be converted to the GMSH file format first. To exchange mesh files between SALOME and GMSH one can use the `*.unv` file format. Finally, the tool GMSH2OGS can be used for the creation of `*.vtu` files (Attention: GMSH2OGS currently accepts the GMSH v2. ASCII format only. It might be also necessary to use the `-e` option to exclude lines before using the mesh files in OGS).

To extract a surface mesh (in order to define appropriate boundary conditions) one can use the tool [ExtractSurface]({{< ref "extract-surface" >}}) for different kinds of 3D meshes. Different tools like [ParaView](https://www.paraview.org/) are suited as well. Care must be taken to make sure that all element types are of dimension `d-1`, where `d` is the dimension of the bulk mesh. Furthermore, one needs to ensure that no degenerated elements are left in the mesh and all nodes are connected appropriately.

Alternatively, the tool [msh2vtu](https://github.com/dominik-kern/msh2vtu) can be used to convert GMSH meshes to VTU and to extract boundary meshes directly from physical groups defined in GMSH.

Heterogeneous fields, e.g., for the use as initial conditions, can be easily generated using the [VTUinterface](https://github.com/joergbuchwald/VTUinterface) Python package.

Finally, use the OGS utility [identifySubdomains]({{< ref "identifySubdomains" >}}) to make sure that the boundary and subdomain meshes are compatible with the bulk mesh in terms of node and element numbering. Using msh2vtu you may skip this step, since all submeshes are extracted from the bulk mesh.
