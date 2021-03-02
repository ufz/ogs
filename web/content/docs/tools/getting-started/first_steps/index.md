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

To set up a model domain, the project file `*.prj` requires additional input files for the mesh and geometry to apply the boundary conditions accordingly.

The mesh must be provided in VTK's `*.vtu` format, whereas the geometry can be provided either in OGS' internal `*.gml` file format (not to be confused with the Geography Markup Language) or as a `*.vtu` file as well, containing boundary elements only (Multiple `*.vtu` files can be provided in the project file using the `<meshes>` tag). We recommend the usage of the first method for simple 2D meshes with constant boundary conditions, whereas more complicated geometries and conditions might require the latter method. One general advantage in the utilization of `*.vtu` files is that they allow a definition of additional field variables at each mesh node/element in order to implement spatially varying boundary conditions in a similar manner as defining inhomogeneous material properties.

The finite-element mesh must be created using external mesh generators. Simple meshes can be created using the [generateStructuredMesh]({{< ref "structured-mesh-generation" >}}) tool. More complicated geometries can be generated and meshed using the [SALOME Platform](https://www.salome-platform.org) or [GMSH](http://gmsh.info).

To use a mesh created by SALOME it must be converted to the GMSH file format first. To exchange mesh files between SALOME and GMSH one can use the `*.unv` file format. Finally, the tool GMSH2OGS can be used for the creation of `*.vtu` files (Attention: GMSH2OGS currently accepts the GMSH v2. ASCII format only. It might be also necessary to use the `-e` option to exclude lines before using the mesh files in OGS).

To extract a surface mesh (in order to define appropriate boundary conditions) one can use the tool [ExtractSurface]({{< ref "extract-surface" >}}) for different kinds of 3D meshes. If you are using a different tool (e.g. paraview), keep in mind that the dimensions of the boundary mesh need to be `d-1`, where `d` is the dimension of the bulk mesh (i.e. `connectivity`, `offsets` and `types` cell properties need to be set correctly in the unstructured grid file).

Finally, the subdomains need to be identified. This can be done by the tool [identifySubdomains]({{< ref "identifySubdomains" >}}) which is also part of the official [OpenGeoSys git repository](https://github.com/ufz/ogs).
