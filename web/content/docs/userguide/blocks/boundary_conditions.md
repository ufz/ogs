+++
date = "2018-02-27T11:00:13+01:00"
title = "Boundary conditions"
author = "Feliks Kiszkurno"
weight = 8
+++
<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

All types of boundary conditions discussed in this article can be defined directly in the project file (more details can be
found [here](/docs/userguide/blocks/process_variables/#boundary-conditions)).
Additionally Dirichlet and Neumann boundary conditions can be defined using [Python boundary conditions functionality](/docs/userguide/features/python_bc/).

## Dirichlet

TODO: add description of Dirichlet boundary condition (including different subvariants) and how it can be used in prj file

### Constrained Dirichlet

### Dirichlet within time interval

### Primary variable constraint Dirichlet boundary condition

### Solution dependent Dirichlet

## Neumann

TODO: add description of Dirichlet boundary condition (including different subvariants) and how it can be used in prj file

Can be used to define [source terms](/docs/userguide/blocks/process_variables/#sources) at the boundary of the mesh.

### Non uniform variable dependent Neumann

### HC non advective free component flow boundary

### Normal traction

## Robin

## Python

Can be used to define [Dirichlet](/docs/userguide/blocks/boundary_conditions/#dirichlet) and [Neumann](/docs/userguide/blocks/boundary_conditions/#neumann) boundary conditions. More details can be found in [this article](/docs/userguide/features/python_bc/).

## Setting source term at a boundary

In some models, the position of the source (of for example heat) can coincide with the boundary condition.
It is possible to use boundary condition to define source.
You need to keep in mind, that in such a case, from the numerical point of view, the source term defined as boundary condition is a boundary condition.
It will be integrated not over the domain (as it would be if it was defined as source term), but over
the boundary condition.
This is why it has to be scaled (see [next section](/docs/userguide/blocks/boundary_conditions/#scaling-a-source)).

Defining it as boundary condition can be done in the same way as it was described in the paragraph above.

### Scaling a source

Neumann boundary condition should be scaled by the surface of the source:

$$ Q_{scaled} = Q_{term} / a $$

where $Q$ is source term, $a$ and resulting unit is power over scaling surface area or length depending on dimension (e.g., $[W / L^{dim-1}]$).
For example if the source is 2D, the scaling factor $a$ will be the length of it's circumference (or the relevant part of it, if
the source is at the boundary of the domain) and similarly in 3D case, the relevant area.
If the mesh is axially symmetric, it has to be considered that it is fully symmetric.
The 2D source will be rotated by $360°$.
For example a rectangular source at the axis of symmetry, will have to be considered as a cylinder.
Then, for scaling the surface of the cylinder has to be used, not the relevant part of the circumference of the rectangle.

If a source is placed on the axially symmetrical boundary and has radius equal to 0, it has to be defined as boundary condition
not as source term.
If radius is non zero, it can be defined using source term.

For more detailed description see [scaling a source term](/docs/userguide/blocks/misc/scaling_source_term/).

TODO: This text is the same as in **Process variables**
