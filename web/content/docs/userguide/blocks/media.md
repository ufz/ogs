+++
date = "2018-02-27T11:00:13+01:00"
title = "Media"
author = "Feliks Kiszkurno"
weight = 4
+++

<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

<!-- TODO: Give links to specific examples. The description found below is still somewhat generic, examples help here. -->

</div>

Inside of this block, all media of the simulation are defined.
There has to be a medium for every value of Material IDs used in the mesh.
Those IDs can be assigned when the mesh is created or converted from the Gmsh
.msh format to the .vtu format via ogstools.

```xml
<media>
    <medium id="MaterialID_matching_mesh_file">

    </medium>
</media>
```

## Phases

A medium can consist of multiple phases.
For example geological porous media would consist of a solid phase (for example gravel, sand, silt, clay, or some mix of the former) and a liquid phase (in most cases water, which fills the pore space either fully or just partially).

In media with multiple phases, each of them can occur at maximum once per medium.
Phases within one medium are distinguished by their type:

```xml
<medium id="0">
    <phases>
        <phase>
            <type>Solid</type>
        </phase>
    </phases>
</medium>
```

If only one phase is considered inside a medium, the definition of that phase can be reduced to placing `<phases/>` as follows:

```xml
<medium id="MaterialID">
    <phases/>
<medium>
```

The available phases depend on the type of the process.

<!-- TODO: Add overview which media are available in different processes -->

## Properties

The tag `<properties> </properties>` either contains all values describing the properties of the medium as a whole or just the properties of a specific phase of the medium.
The properties, which describe the medium as a whole, have to be placed on the same level as phases:

```xml
<medium id="0">
    <phases>
        ...
    </phases>
    <properties>
        Properties describing the whole medium
    </properties>
</medium>
```

The properties belonging to a specific phase within a medium, have to be placed inside `<phase> </phase>` tag:

```xml
<medium id="0">
    <phases>
        <phase>
            <type>Solid</type>
            <properties>
                Properties describing the solid phase medium 0
            </properties>
        </phase>
    </phases>
</medium>
```

The properties share very similar structure with [parameters](/docs/userguide/blocks/parameters), but they are not interchangeable.
The parameters can represent time and space dependent values, whereas the medium properties are representing material laws.
The latter can refer to the parameters and therefore the parameters are more basic in their nature.

Following basic types of properties are available:

- [Constant](/docs/userguide/blocks/media/#constant)
- [Linear](/docs/userguide/blocks/media/#linear)
- [Function](/docs/userguide/blocks/media/#function)
- [Curve](/docs/userguide/blocks/media/#curve)

They can be used by the user to define the properties in a way, that is specific to their simulation.

There are more general properties available.
They are described in section [Other types of properties](/docs/userguide/blocks/media/#other-types-of-properties).

The types linear, function, and curve can depend on a set of variables listed in [MPL->VariableType.h](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/MaterialLib/MPL/VariableType.h).
Keep in mind that not all of those variables will be available in all the processes.

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

For vectorial quantities, please list the individual entries.
The order of the entries has to follow the order of the x, y and z components of the vector.
`<values>0.5 1 2</values>`, will be interpreted as a 3D vector with components for x, y and z.

### Linear

The type linear can be used to declare parameters that vary linearly depending on the value of another parameter.

Apart from the standard tags, which are required: `<name> </name>` and `<type> </type>`, `<reference_value> </reference_value>`
has to be provided, too.

Linearly varying parameters can be a function of space (x, y, and z) and time (t).
Linear dependency on each independent variable requires `<reference_condition> </reference_condition>` and `<slope> </slope>`.

Linear property is defined as follows:

$$
y(x_{i}) = y_{\textrm{ref}} \left(1 + \sum_{i=1}^{n} m_{i} (x_{i} - x_{i,\textrm{ref}})\right)
$$

where

- $x_{i}$ can be a number of dependent variables, for instance temperature, pressure
- $y_{\textrm{ref}}$ is a reference value, for instance reference density
- $m_{i}$ is a value influencing the slope of the linear relationship with respect to dependent variable
- $x_{i, \textrm{ref}}$ is a reference condition with respect to dependent variable, for instance reference temperature, reference pressure

Defining it as type `linear` would look like this in project file:

```xml
<property>
    <name>y</name>
    <type>Linear</type>
    <reference_value>y_ref</reference_value>
    <independent_variable>
        <variable_name>dependent_variable_x</variable_name>
        <reference_condition>x_ref</reference_condition>
        <slope>m</slope>
    </independent_variable>
</property>
```

where the reference value $y_{ref}$ is defined as follows:

$$
y_{ref}=y(x_{ref})
$$

For this block the equation would simplify to:

$$
y(x_{i}) = y_{\textrm{ref}} \left(1 + m (x - x_{\textrm{ref}})\right)
$$

as there is only one dependent variable.

A more realistic example can be found in benchmark [A2](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/Tests/Data/HydroMechanics/A2/A2.prj), where the fluid density linearly but independently depends on pressure and temperature:

```xml
<property>
    <name>density</name>
    <type>Linear</type>
    <reference_value>1200</reference_value>
    <independent_variable>
        <variable_name>temperature</variable_name>
        <reference_condition>298.15</reference_condition>
        <slope>-6.0e-4</slope>
    </independent_variable>
    <independent_variable>
        <variable_name>phase_pressure</variable_name>
        <reference_condition>4e6</reference_condition>
        <slope>0.5e-9</slope>
    </independent_variable>
</property>
```

The slope value provided in the `<independent_variable> </independent_variable>` block for temperature is the thermal expansion coefficient.

The definition of density provided in the snippet above can be expressed by following equations:
$$
\rho_{0} (T_{0}=298.15, p_{0}=4\cdot 10^{6})=1200
$$

$$
\rho (T) = \rho_{0} (1 + (-6 \cdot 10^{-4})\cdot (T-T_{0}))
$$

$$
\rho (p) = \rho_{0} (1 + (0.5 \cdot 10^{-9}) \cdot (p-p_{0}))
$$

### Function

The type `function` extends the basic structure of [constant type](/docs/userguide/blocks/parameters/#constant).
There are two major differences syntax-wise though.
The expression that defines a function has to be placed inside of the `<expression> </expression>` tag nested inside of the `<value> </value>` tag.
Additionally, derivatives with respect to applicable process variables have to be provided inside a template as follows:

```xml
<dvalue>
    <variable_name>primary variable derived over</variable_name>
    <expression>derivative of function w.r.t. specific primary variable </expression>
</dvalue>
```

#### Convention & Syntax

Powers are indicated by "^", e.g.: 2^3=8.

However, following standard conventions, powers to the basis of 10 may also be expressed as  `0.002=2e-3=2e-03`.
All three expression can be used in OpenGeoSys and are obviously equivalent.

C++ functions (for example std::pow) cannot be called from the expression tag.

The [`exprtk`](http://www.partow.net/programming/exprtk/) library interpreting such expressions is used in OGS and the full list of available functions and control structures can be found on its web page.

#### Example

Consider the following function defining viscosity of water depending on temperature:

$$ \mu = ( -0.0002 \cdot T^3 + 0.05 \cdot T^2 - 4 \cdot T + 178) \cdot 10^{-5}$$

and its partial derivative with respect to temperature:

$$ \frac{\partial \mu}{\partial T} = ( -0.0006 \cdot T^2 + 0.1 \cdot T - 4) \cdot 10^{-5} $$

With the two equations above, the full implementation of water viscosity depending on phase pressure and temperature will look
as follows:

```xml
<property>
    <name>viscosity</name>
    <type>Function</type>
    <value>
        <expression>(-0.0002*(temperature-273.15)^3+0.05*(temperature-273.15)^2-4*(temperature-273.15)+178)*1e-5</expression>
    </value>
    <dvalue>
        <variable_name>temperature</variable_name>
        <expression>(-0.0006*(temperature-273.15)^2+0.1*(temperature-273.15)-4)*1e-5</expression>
    </dvalue>
</property>
```

There are limitations of what variables can be used inside of the `<expression> </expression>` tags.
Only the ones related to the process variables can be called. The values defined for example in the `parameter` block are out of reach.

Notice, that OpenGeoSys assumes self consistent set of units.
In the code snippet above temperature has to be converted from kelvins (used by other temperature dependent parameters in the `.prj`-file from which this snippet comes) to degrees Celsius required by the provided formula.
OpenGeoSys doesn't have an internal unit management system.
It is your task to ensure, that the units used are consistent!

### Curve

Curve is a special type, as it requires the presence of another block in the project file.
Whether the parameter is defined in the Parameters block, or anywhere else in the project file, type curve will always refer to the curve defined in [curve block](/docs/userguide/blocks//curves/).
Regardless, of whether there is only one curve defined or multiple, they are always identified by the content of the tag `<name> </name>`.

![Relation between parameter of type curve and curves block](/docs/userguide/blocks/figures/curve.svg)

Parameters of the type `curve` can be defined using the following template:

```xml
<property>
    <name>property_name</name>
    <type>Curve</type>
    <independent_variable>curve_name_coords</independent_variable>
    <curve>curve_name</curve>
</property>
```

where independent variable refers the values in tag `<coords> </coords>` in block `<curves> </curves>`.
For more details see [example of curve definition](/docs/userguide/blocks/curves/#example).

Type "Curve" is different to the one "Properties" tag.
In the "Expression" tag only numerical values and x, y, z, t variables and curves can be used.

## Other types of properties

Some properties being an outcome of other models are implemented in OpenGeoSys and may be used as well.
In this section some examples will be discussed, however, more models than discussed are available in OGS.

### Water properties - IAPWS97 model

The water property model [IAPWS97](http://www.iapws.org/relguide/IF97-Rev.html) ([International Association for the Properties of Water and Steam](http://www.iapws.org)) is implemented in OpenGeoSys.
Various properties of water can be called by providing an appropriate tag:

- density `<type>WaterDensityIAPWSIF97Region1</type>`
- viscosity `<type>WaterViscosityIAPWS</type>`
- thermal conductivity `<type>WaterThermalConductivityIAPWS</type>`

Those tags can be used only on the phase level.
Please see the following example:

```xml
<phase>
    <type>AqueousLiquid</type>
    <properties>
        <property>
            <name>density</name>
            <type>WaterDensityIAPWSIF97Region1</type>
        </property>
        <property>
            <name>viscosity</name>
            <type>WaterViscosityIAPWS</type>
        </property>
        <property>
            <name>viscosity</name>
            <type>WaterThermalConductivityIAPWS</type>
        </property>
    </properties>
</phase>
```

For those properties SI units are used, remember that OpenGeoSys doesn't convert units internally.

### Supported types of properties and special cases

<!-- TODO: Add content -->
