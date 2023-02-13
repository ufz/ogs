+++
date = "2018-02-27T11:00:13+01:00"
title = "Parameters"
author = "Feliks Kiszkurno"
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

The possible content of this block in not limited to physical properties of materials used in the experiment but also to values of [boundary](/docs/userguide/blocks/process_variables/#boundary-conditions) and [initial conditions](/docs/userguide/blocks/process_variables/#initial-conditions) or [source terms](/docs/userguide/blocks/process_variables/#sources), physical constants, etc.

## Where can the parameters be used?

The parameters defined in this block can be used in blocks:

- [Media](/docs/userguide/blocks/media/)
- [Process variables](/docs/userguide/blocks/process_variables/)

## Parameters vs properties

TODO: describe differences in access to parameters and properties

## How to define a parameter?

To create a parameter within the `<parameters> </parameters>` tag, following template can be used:

```xml
<parameter>
    <name>...</name>
    <type>...</type>
</parameter>
```

Tags `<name> </name>` and `<type> </type>` are mandatory and define human-readable name of specific parameter and declare one of available types with which it will be defined.
Other tags depend on what is the content of `<type> </type>`.
There are following types available:

- [Constant](/docs/userguide/blocks/parameters/#constant)
- [CurveScaled](/docs/userguide/blocks/parameters/#curvescaled)

Each of them will be discussed below.
The same types can be used to define media properties in the [media block](/docs/userguide/blocks/media/).

### Constant

This is the most basic type.
It only requires `<value> </value>` tag additionally where the value of the parameter is provided as a number.
It will not change throughout the experiment.
For example:

```xml
<parameter>
    <name>earth_acceleration</name>
    <type>Constant</type>
    <value>9.81<value>
</parameter>
```

For vectorial quantities use plural, `<values>0.5 1 2</values>`, for example will be interpreted as a 3d vector:

```xml
<parameter>
    <name>some_anisotropic_parameter</name>
    <type>Constant</type>
    <value>0.5 1 2</value>
</parameter>
```

### CurveScaled

It requires tags: `<curve> </curve>` and `<parameter> </parameter>`.
The first one contains name of a curve defined in the `<curves> </curves>` block (for more detail see [curves](/docs/userguide/blocks/curves/)).
It is not the same type as [Curve](/docs/userguide/blocks/media/#curve) discussed in Media block.

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
