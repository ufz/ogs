+++
date = "2018-02-27T11:00:13+01:00"
title = "Meshes"
author = "Feliks Kiszkurno"
weight = 2
+++

Here, a minimum of two files has to be provided.
One for the mesh of the domain and one for the boundaries. However, arbitrary many mesh files can be provided.
Generally, if the boundary conditions vary between boundaries, there should be one mesh file provided for each set of boundaries
sharing the same conditions. Alternatively, boundary conditions can be manipulated via Python, too.

<div class = note>

TODO: Give examples for the both files named above and how to incorporate them. Give a link to a tutorial, how to build mesh-files as well as boudary-files.

</div>

```xml
  <mesh>mesh_file_name.mesh_extension</mesh>
```

## Axial symmetry

Additionally it can be specified, that the provided mesh should be axially symmetric:

```xml
  <mesh axially_symmetric="true">
      mesh_file_name.mesh_extension
  </mesh>
```

If this flag is set to "True", the mesh will be "fully", i.e., 360° symmetric.
This also has an influence on the appropriate scaling of the source term (see [Source term](/docs/userguide/blocks/processes/#source-term)).

OGS follows a standard axis convention: XY is a horizontal plane and Z is vertical.
By default the rotation will be around the Y-axis. Moreover, see [specific body force](/docs/userguide/blocks/processes/#specific-body-force) for details on including the gravitational force in accordance with the axis definitions.
