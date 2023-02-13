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

</div>

Inside of this block, all media in the simulation are defined.
 There has to be a medium for every value of Material IDs used in the mesh.
 Those IDs are assigned when mesh is created or converted by msh2vtu script.

```xml
<media>
    <medium id="MaterialID_matching_mesh_file">

    </medium>
</media>
```

## Phases

Medium can consist of multiple phases.
For example porous geological would consist of solid phase (for example clay) and liquid phase (e.g.: water in the pores).

In media with multiple phases, each of them can occur maximal once per medium.
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

TODO: Add overview of what mediums are available in different processes

## Properties

The tag `<properties> </properties>` contains all values describing the properties of the medium as whole or of the specific phase of the medium. 
The properties describing medium as a whole have to be placed on the same level as phases:

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

The properties belonging to a phase within a medium, have to be placed inside `<phase> </phase>` tag:

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

The properties share very similar structure with [parameters](/docs/userguide/blocks/parameters/), but they are not interchangeable. 

Following basic types of properties are available:

- [Constant](/docs/userguide/blocks/media/#constant)
- [Linear](/docs/userguide/blocks/media/#linear)
- [Function](/docs/userguide/blocks/media/#function)
- [Curve](/docs/userguide/blocks/media/#curve)

They can be used by the user to define the properties in a way, that is specific to their simulation.

There are more general properties available.
They are described in section [Other types of properties](/docs/userguide/blocks/media/#other-types-of-properties).

Generally, it is the most safe to use "Constant" type for properties and when it is not sufficient, type "Parameter" can be used as fallback.
Still, there are some limitations to what types of parameter can be used in different processes.

Differently than the parameters in the Parameter block, in the Media block parameters can depend on more variables than x, y, z and t.

Type linear, function and curve can depend on a set of variable listed in [MPL->VariableType.h](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/MaterialLib/MPL/VariableType.h):

- capillary_pressure
- concentration
- density
- displacement
- effective_pore_pressure
- enthalpy_of_evaporation
- equivalent_plastic_strain
- grain_compressibility
- liquid_phase_pressure
- liquid_saturation
- mechanical_strain
- molar_mass
- molar_mass_derivative
- molar_fraction
- phase_pressure
- porosity
- solid_grain_pressure
- stress
- temperature
- total_strain
- total_stress
- transport_porosity
- vapour_pressure
- volumetric_strain

Keep in mind that not all of those variables will be available in all the processes. For example, in THM there is phase_pressure, but not liquid_phase_pressure.

### Constant

This is the most basic type. It only requires `<value> </value>` tag additionally where the value of the parameter is provided as a number. It will not change throughout the experiment. For example:

```xml
<parameter>
    <name>earth_acceleration</name>
    <type>Constant</type>
    <value>9.81<value>
</parameter>
```

For vectorial quantities use plural, `<values>0.5 1 2</values>`, for example will be interpreted as a 3d vector.

### Linear

Type linear, can be used to declare parameters that vary linearly depending on value of another parameter.

Apart from standard required tags: `<name> </name>` and `<type> </type>`, `<reference_value> </reference_value>` has to be provided.

Linear parameter can depend on values describing position (x, y and z) and time (t).
Linear dependency on each independent variable requires `<reference_condition> </reference_condition>` and `<slope> </slope>`.

Consider simple linear equation:
$$
y(x)=ax+b
$$
let's say that there is a parameter $y$ that depends linearly on time:
$$
y(T)=a*t+y_0
$$
Defining it as linear type would look like this in project file:

```xml
<property>
    <name>y</name>
    <type>Linear</type>
    <reference_value>y_0</reference_value>
    <independent_variable>
        <variable_name>time</variable_name>
        <reference_condition>t_0</reference_condition>
        <slope>a</slope>
    </independent_variable>
</property>
```

where reference value $y_0$ is defined as follows:
$$
y_0=y(t_0)
$$

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

The definition of density provided in the snippet above can be expressed by following equations:
$$
\rho(T=298.15, p=4e6)=1200
$$
$$
\rho(T)=(-6*10^{-4})*T+1200
$$
$$
\rho(p)=(0.5*10^{9})*p+1200
$$

### Function

Type function extends the basic structure of [constant type](/docs/userguide/blocks/parameters/#constant).
There are two major differences though.
The expression that defines a function has to be placed inside of the `<expression> </expression>` tag nested inside of the `<value> </value>` tag.
Additionally, derivations with respect to applicable process variables have to be provided inside following this template:

```xml
<dvalue>
    <variable_name>primary variable derived over</variable_name>  
    <expression>expression of derivation of function w.r.t. specific primary variable </expression>
</dvalue>
```

#### Convention & Syntax

Powers are indicated by "^", e.g.: 2^3=8.

In scientific notation, power can be appended by leading zeros. So `0.002=2e-3=2e-03`, those three expression will evaluate to the same value in OpenGeoSys.

C++ functions (for example std::pow) cannot be called from the expression tag.

The [`exprtk`](http://www.partow.net/programming/exprtk/) library interpreting such expressions is used in OGS and the full list of available functions and control structures can be found on its web page.

#### Example

Consider function defining viscosity of water depending on temperature:

$$ \mu = ( -0.0002 * T^3 + 0.05 * T^2 - 4 * T + 178) * 10^{-5}$$

and its partial derivative with respect to temperature:
$$ \frac{\partial \mu}{\partial T} = ( -0.0006 * T^2 + 0.1 * T - 4) * 10^{-5} $$

With the two equations above, the full implementation of water viscosity depending on phase pressure and temperature will look as follows:

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

There are limitations of what variables can be used inside of the `<expression> </expression>` tags. Only the ones related to the process variables can be called. The values defined for example in parameter block are out of reach.
TODO: Check if this is correct

For example if ```thermal_conductivity``` is defined and density is provided as Function, the dvalue of density with respect to temperature will be ignored.

Notice, that OpenGeoSys assumes self consistent set of units.
In the code snippet above temperature has to be converted from kelvins (used by other temperature dependent parameters in the project from which this snippet comes) to degrees Celsius required by the provided formula.
OpenGeoSys doesn't have an internal unit management system.
It is your task to ensure, that the units used are consistent and compatible!

### Curve

Curve is a bit special type, as it requires presence of another block in the project file.
Whether the parameter is defined in Parameters block, or anywhere else in the project file, type curve will always refer to curve defined in [curve block](/docs/userguide/blocks//curves/).
Regardless, of whether there is only one curve defined or multiple, they are always identified by the content of tag `<name> </name>`.

![Relation between parameter of type curve and curves block](/docs/userguide/blocks/figures/curve.svg)

Parameter of type curve can be defined using the following template:

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

Some properties that are valid for different simulation and can be reused are implemented in OpenGeoSys.
In this section they will be discussed, but this list is not complete. 

### Water properties - IAPWS97 model

The water property model [IAPWS97](http://www.iapws.org/relguide/IF97-Rev.html) ([International Association for the Properties of Water and Steam](http://www.iapws.org)) is implemented in the OpenGeoSys.
Various properties of water can be called by providing an appropriate tag:

- density `<type>WaterDensityIAPWSIF97Region1</type>`
- viscosity `<type>WaterViscosityIAPWS</type>`
- thermal conductivity `<type>WaterThermalConductivityIAPWS</type>`

Those tags can be used only on the phase level. Please see following example:

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

For those properties SI units are used, OpenGeoSys doesn't convert units internally.
Users need to ensure, that the parameters that are provided in the project file use consistent and correct units.

### Supported types of properties and special cases

TODO: Add content
