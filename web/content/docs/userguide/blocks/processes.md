+++
date = "2018-02-27T11:00:13+01:00"
title = "Processes"
author = "Feliks Kiszkurno"
weight = 3
+++
<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

Process block is defined with `<process> </process>` tags.

This section defines the process, that is modelled in the project.
It has to contain 'name' and 'type' tags:

```xml
<name>THM</name>
<type>THERMO_HYDRO_MECHANICS</type>
```

Name can be defined freely, type has to contain one of following strings [list of available processes](https://doxygen.opengeosys.org/d5/d98/ogs_file_param__prj__processes__process.html). More details about some of them can be found in the [Process-dependent configuration](/docs/processes/heat-transport/heat_transport_bhe/) section.

## Integration order

Tag `<integration_order> </integration_order>` is mandatory and defines the integration order of the Gauss-Legendre integration over an element.

## Process variables

An important part of this section is defined in "process\_variables" tag. It is important, because in later parts of the project files (e.g.: in definition of errors). Those variables and their order are specific to each type of process. The order in which those variables have been defined in each process plays an important role.
For more details see [time loop](/docs/userguide/blocks/time_loop/) and [linear solvers](/docs/userguide/blocks/linear_solvers/).
Information of which variables are required for which process can be found in appropriate sections of the Doxygen documentation.
For example, for THM process following parameters are required:

```xml
<process_variables>
    <temperature>temperature</temperature>
    <pressure>pressure</pressure>
    <displacement>displacement</displacement>
</process_variables>
```

## Secondary variables

This tag is optional.
It allows to define which secondary variables will be computed and saved in the output files.
Which secondary variables are available is dependent on a specific process.
They can be defined as follows:

```xml
<secondary_variables>
    <secondary_variable internal_name="internal_name" output_name="output_name"/>
</secondary_variables>
```

where internal name has to match the name of one of available secondary variables and output name can be defined by the user.
Output name will be used as the name of the field into which specific variable will be written into in output files.

## Error tolerances

In this section, within `<convergence\_criterion>` tag, relative tolerances for error have to be defined - with tag `<reltols>`.
Each value in this tag defines the tolerance for errors with respect to one process variables.
The order of variables is process specific.
For the process variables listed above, relative error tolerances can be defined as follows:

```xml
<reltols>
    reltol_temp
    reltol_press
    reltol_disp_0
    reltol_disp_1
</reltols>
```

The order can differ based on the order in which they are defined in a process and on dimensionality of the process (*e.g.* number of components in displacement, or number of chemical constituents).
Keep in mind, that some process variables have more than one values, as displacement in the example above.
In such a case, a matching number of reltols has to be defined.

## Constitutive relations

Constitutive relation can be one of the [existing relations](/docs/userguide/blocks/misc/constitutive_relations/) implemented in OpenGeoSys or it can be defined by user using [MFront](/docs/userguide/features/mfront/).
They are used with one of the following mechanical processes:

* Hydro Mechanics
* Phase Field
* Richards Mechanics
* Small Deformation
* Thermo Mechanics
* Thermo Hydro Mechanics
* Thermo Richards Mechanics
* TH2M

To define constitutive relation, tags `<constitutive_relation> </constitutive_relation>` are used.
The fixed, minimum requirement is presence of `<type> </type>` tag. Other tags depend on the chosen relation.

There can be more than one constitutive relation used in one project file but only one for each material id, similar to how each medium is defined in [media block](/docs/userguide/blocks/media/).
OpenGeoSys selects constitutive relation for a material id based on the content of the attribute id in tags `<constitutive_relation> </constitutive_relation>` and `<medium> </medium>`.

## Source terms

There are four available types of source terms:

* Line
* Nodal
* Volumetric
* Python

 The source term has to be scaled by the volume.  

## Specific body force

 The tag `<specific_body_force>` can be used to define gravitational force. 
 This tag is optional.
 By default it is not considered in the simulation, it has to be be explicitly enabled by setting an appropriate value in project file.
 Presence of this tag, should be reflected in initial stress conditions.
 It can improve the convergence in the first time step, if the correct stress field is defined.

```xml
<process>
    ...
    <specific_body_force>in_X_direction in_Y_direction in_Z_direction
    </specific_body_force>
    ...
</process>
```

In 3D case it has to be defined on Z-axis as negative value: $-9.81$.
So the code in project file would look like this:

```xml
<process>
    ...
    <specific_body_force>0 0 -9.81</specific_body_force>
    ...
</process>
```

**Warning!** If specific body force is added, it can be reflected with the initial state of stress.

Stress is a a bit unique, as initial condition for stress are not defined in block "Process variables", but in "Process".
It requires setting up a constant or function definition in "Parameters" block where initial stress is defined.
It has to match the dimension of experiment.
OGS expects it to be effective stress (total stress with pore pressure subtracted), this is in contrast to total stress that is expected to be provided in boundary condition.
If meshes are defined as axially symmetrical, the stress has to be provided in polar or cylindrical coordinates.

## Jacobian Assembler

TODO: Explanations for each type.

The global non-linear equation system can be solved either with Picard fix-point iterations or a Newton scheme. (TODO: Reference NLS scheme)
For the latter the process has to provide an analytical Jacobian.
If the analytical Jacobian was not implemented, one can use a quasi-Newton scheme where the Jacobian is approximated by a central differences or forward differences scheme.

For the Picard fix-point iterations this section is not used.
