+++
date = "2018-02-27T11:00:13+01:00"
title = "Process variables"
author = "Feliks Kiszkurno"
weight = 7
+++

All "process parameters" defined in the "Process" section have to be addressed here (see [processes](/docs/userguide/blocks/processes/)).

```xml
<process_variable>
    <name>"process_variable"</name>
    <components>integer</components>
    <order>integer</order>
    <initial_condition>variable_name</initial_condition>
    <boundary_conditions>
        <boundary_condition>
            <type>Dirichlet</type>
            <mesh>Temp_Dis_mesh_boundary</mesh>
            <component>1</component>
            <parameter>dirichlet0</parameter>
        </boundary_condition>
    </boundary_conditions>
</process_variable>
```

The block `<process_variable>` has be defined separately for each process variable defined in [Block "Process/Process variables"](/docs/userguide/blocks/processes/#process-variables).

TODO: Explain the little example above in a little bit more detail or give a link to a tutorial.

There are three major blocks in this one.
Two of them are mandatory, in following order:

- [Initial conditions](/docs/userguide/blocks/process_variables/#initial-conditions)
- [Boundary conditions](/docs/userguide/blocks/process_variables/#boundary-conditions)

and optionally:

- [Source terms](/docs/userguide/blocks/process_variables/#sources)

can be defined.

The tag <components> refers to how many directional components a variable has, for example, displacement in a 2D THM process
will have 2 components (for x and y directions).

`<initial_conditions>`tag should contain the name of variable defined in [parameters block](/docs/userguide/blocks/parameters/).

## Initial conditions

Initial conditions define the state of the primary variables in the domain at $t=0$.
Each process will require a different set of them, see [process variables](/docs/userguide/blocks/processes/#process-variables)
for more details.

For example, the value of initial temperature can be defined in [parameters block](/docs/userguide/blocks/parameters/) and
called from `<process_variable> </process_variable>` by the content of `<name> </name>` tag.
See following example:

```xml
...
    <name>temperature</name>
    <initial_condition>initial_temperature</initial_condition>
...
...
<parameter>
    <name>initial_temperature</name>
    <type>Constant</type>
    <value>273.15</value>
</parameter>
```

## Boundary conditions

The two most commonly used boundary condition types are [Dirichlet](/docs/userguide/blocks/boundary_conditions/#dirichlet) and [Neumann](/docs/userguide/blocks/boundary_conditions/#neumann).
The first one can be used to set boundary conditions to a set value, Neumann can be used for flow (e.g., water or energy).
For more details on both issues, see [sources](/docs/userguide/blocks/process_variables/#sources).
For simulations which require more complicated values at boundaries (for example spatially variable), [Python boundary conditions](/docs/userguide/blocks/boundary_conditions/#python) can be used.

It is possible to define a source at the boundaries of the simulation domain using boundary conditions.
For more details, see [Setting source term at a boundary](/docs/userguide/blocks/boundary_conditions/#setting-source-term-at-a-boundary)

## Sources

Source terms are optional.
There are two ways how sources can be defined in the project file: as boundary condition or as source term.

### Source term

There are following types available: nodal, volumetric.

```xml
<source_term>
    <geometrical_set>geometry</geometrical_set>
    <geometry>inner</geometry>
    <type>Nodal</type>
    <parameter>pressure_source_term</parameter>
</source_term>
```

### Source as boundary condition

Defining it as boundary condition can be done in the same way as it was described in the paragraph above.

### Scaling a source

Neumann boundary condition should be scaled by the surface of the source:

$$ Q_{scaled} = Q_{term} / a $$

where $Q$ is source term, $a$ and resulting unit is power over scaling surface area or length depending on dimension (e.g., $[W / L^{dim-1}]$).
For example if the source is 2D, the scaling factor $a$ will be the length of it's circumference (or the relevant part of it,
if the source is at the boundary of the domain) and similarly in 3D case, the relevant area.
If the mesh is axially symmetric, it has to be considered that it is fully symmetric.
The 2D source will be rotated by $360°$.
For example a rectangular source at the axis of symmetry, will have to be considered as a cylinder.
Then, for scaling the surface of the cylinder has to be used, not the relevant part of the circumference of the rectangle.

If a source is placed on the axially symmetrical boundary and has radius equal to 0, it has to be defined as boundary condition
not as source term.
If radius is non zero, it can be defined using source term.

For more detailed description see [scaling a source term](/docs/userguide/blocks/misc/scaling_source_term/).
