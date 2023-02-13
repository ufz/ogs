+++
date = "2018-02-27T11:00:13+01:00"
title = "Meshes"
author = "Feliks Kiszkurno"
weight = 2
+++

Minimum of two files has to be provided.
One for the mesh of the domain and one for the boundaries. Arbitrary many mesh files can be provided.
Generally, if the boundary conditions vary between boundaries, there should be one mesh file provided for each set of boundaries sharing the same conditions.
Alternatively Python boundary conditions can be used.

```xml
  <mesh>mesh_file_name.mesh_extension</mesh>
```

## Axial symmetry

Additionally it can be specified, that provided mesh is axially symmetric:

```xml
  <mesh axially_symmetric="true">
      mesh_file_name.mesh_extension
  </mesh>
```

if it is set to "True", the mesh will be "fully" - 360° symmetric.
This is relevant for scaling of the source term (see [Source term](/docs/userguide/blocks/processes/#source-term)).

OGS follows standard axis convention: XY is a horizontal plane and Z is vertical.
By default the rotation will be around Y-axis. See [specific body force](/docs/userguide/blocks/processes/#specific-body-force) for details on including gravitational force.
