+++
date = "2025-08-29T09:00:13+01:00"
title = "Parameters"
author = "Feliks Kiszkurno, Florian Berger"
weight = 5
+++

<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

This block contains various parameters which are used by other blocks within the project file.

The possible content of this block in not limited to physical properties of materials used in the experiment but also to values
of [boundary](/docs/userguide/blocks/process_variables/#boundary-conditions) and [initial conditions](/docs/userguide/blocks/process_variables/#initial-conditions) or [source terms](/docs/userguide/blocks/process_variables/#sources), physical constants, etc.

## Where can the parameters be used?

The parameters defined in this block can be used in blocks:

- [Media](/docs/userguide/blocks/media/)
- [Process variables](/docs/userguide/blocks/process_variables/)

## Parameters vs properties

<!-- TODO: describe differences in access to parameters and properties -->

## How to define a parameter?

To create a parameter within the `<parameters> </parameters>` tag, following template can be used:

```xml
<parameter>
    <name>...</name>
    <type>...</type>
</parameter>
```

Tags `<name> </name>` and `<type> </type>` are mandatory and define a human-readable name of a specific parameter and declare
one of the available types with which it will be defined.
Other tags depend on what is the content of `<type> </type>`.
There are the following types available:

- [Constant](/docs/userguide/blocks/parameters/#constant)
- [CurveScaled](/docs/userguide/blocks/parameters/#curvescaled)
- [Function](/docs/userguide/blocks/parameters/#Function)
- [MeshElement](/docs/userguide/blocks/parameters/#meshelement)
- [MeshNode](/docs/userguide/blocks/parameters/#meshelement)
- [RandomFieldMeshElement](/docs/userguide/blocks/parameters/#RandomFieldMeshElement)
- [TimeDependentHeterogeneous](/docs/userguide/blocks/parameters/#TimeDependentHeterogeneous)

Each of them will be discussed below.
The same type can be used to define media properties in the [media block](/docs/userguide/blocks/media/).

### Constant

This is the most basic type.
It only requires the `<value> </value>` tag additionally where the value of the parameter is provided as a number.
It will not change throughout the experiment.
For example:

```xml
<parameter>
    <name>earth_acceleration</name>
    <type>Constant</type>
    <value>9.81<value>
</parameter>
```

For vectorial quantities, the individual entries need to be listed, e.g., `<values>0.5 1 2</values>`.
This example will be interpreted as a 3D vector.

```xml
<parameter>
    <name>some_anisotropic_parameter</name>
    <type>Constant</type>
    <value>0.5 1 2</value>
</parameter>
```

<!-- TODO: This is already (partially) described in the section **Media** -->

### CurveScaled

It requires the tags `<curve> </curve>` and `<parameter> </parameter>`.
The first one contains the name of a curve defined in the `<curves> </curves>` block (for more detail see [curves](/docs/userguide/blocks/curves/)).
Note that this is not the exact same type as [Curve](/docs/userguide/blocks/media/#curve) discussed in Media block.

Following variables are accessible from CurveScaled:

- spatial: x, y, z
- temporal: t

```xml
<parameter>
    <name>p_Dirichlet_right</name>
    <type>CurveScaled</type>
    <curve>p_Dirichlet_right_temporal</curve>
    <parameter>t</parameter>
</parameter>
```

```xml
<curve>
    <name>p_Dirichlet_right_temporal</name>
    <coords>0.0 10.0</coords>
    <values>0.0 10.0</values>
</curve>
```

### Function

It is also possible to define parameters using a function, enabling spatially varying or time-dependent behavior of the parameter.
These functions are defined in the necessary `<expression></expression>` subtag.
Valid function variables are `x, y, z and t`, so the spatial variable and the time.
You can even incorporate predefined curves from the [curves section](/docs/userguide/blocks/curves/) of the project file.

```xml
<parameter>
    <name>density</name>
    <type>Function</type>
    <expression>2000 - CurveName(x^2 + y^2) * TimeDependentCurve(t)</expression>
</parameter>
```

### Group

The group type parameter is defined by a specific value or values for each of the group ids.
The `<group_id_property></group_id_property>` tag specifies the name of the data array in the mesh that define the group.
In the `<index_values></index_values>` subtags, the values of the parameters in those group ids are defined.
Inside this tag you will find the `<index></index>` subtag, that represents the values of the group id, while the `<values></values>` or `<value></value>` subtags define the values of the parameter for fields where aforementioned group id matches the index.

```xml
<parameter>
    <name>K</name>
    <type>Group</type>
    <group_id_property>MaterialIDs</group_id_property>
    <index_values>
        <index>0</index>
        <values>1 1 1</values>
    </index_values>
    <index_values>
        <index>1</index>
        <values>1 1 0.1</values>
    </index_values>
    <index_values>
        <index>2</index>
        <values>0.1 0.1 1</values>
    </index_values>
    <use_local_coordinate_system>true</use_local_coordinate_system>
</parameter>
```

<h3 id = "meshelement"> MeshElement and MeshNode</h3>

Parameters of type "MeshElement" or "MeshNode" are defined by a data array in a given mesh.
The difference between "MeshElement" and "MeshNode" type is the underlying data type.
While "MeshElements" use cell labels, "MeshNode" uses point labels.
Which to use is dependent on the parameter or the provided data.
Usually primary variables use node data while cell data is used for secondary variables.
See [process variables](/docs/userguide/blocks/processes/#process-variables) for more information on variable types.

The mesh must be loaded in the [mesh section](docs/userguide/blocks/Meshes) of the project file.
To specify the data array which contains the desired information, use the `<field_name></field_name>` subtag.
This field must be present in the mesh defined in the `<mesh></mesh>` subtag.
If the`<mesh></mesh>` tag is not set, the first defined mesh will be used.

```xml
<parameter>
    <name>K_rho_over_mu__eff</name>
    <type>MeshElement</type>
    <field_name>K_rho_over_mu__eff</field_name>
</parameter>
<parameter>
    <name>mass_flux</name>
    <type>MeshNode</type>
    <mesh>inhomogeneous_permeability_top</mesh>
    <field_name>mass_flux</field_name>
</parameter>
```

### RandomFieldMeshElement

Defines a parameter with random values in each mesh element with a uniform distribution in a given range.
The range can be defined in the `<range></range>` subtag.
A given seed (`<seed></seed>` subtag) for initializing the random number generator allows for repeatability of the random field.
The `<field_name></field_name>` tag defines the name of the output variable that will contain this parameter.

```xml
<parameter>
    <name>phi</name>
    <type>RandomFieldMeshElement</type>
    <field_name>phi_xy</field_name>
    <range>0 3.1415926535</range>
    <seed>20210422</seed>
</parameter>
```

### TimeDependentHeterogeneous

<div class="note">

Work in progress.
This section is not yet documented.

</div>
