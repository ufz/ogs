+++
date = "2018-02-27T11:00:13+01:00"
title = "Constitutive relations"
author = "Feliks Kiszkurno"
weight = 14
+++
<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

## Overview

## CreepBGRa
TODO: Add content

## Ehlers
TODO: Add content

## Linear Elastic Isotropic
TODO: Add content

## Linear Elastic Orthotropic
TODO: Add content

## Lubby2
TODO: Add content

## MFront

This section only describes how to use an MFront constitutive relation inside of the project file. For details on how to create such a constitutive relation, please see [MFront](/docs/userguide/features/mfront/) section.

In the project file, MFront model has to be selected as constitutive model in `<process> </process>` tag. See following example:

```xml
<constitutive_relation id="0">
    <type>MFront</type>
    <behaviour>ModelName</behaviour>
    <material_properties>
        <material_property name="MatPropName" parameter="ParameterName"/>
    </material_properties>
</constitutive_relation>
```

In this example ModelName has to match one of the MFront files in the OpenGeoSys source code. [Attribute](/docs/userguide/basics/project_file_intro/#attributes) "name" defines name of one of the material properties in the MFront file (see [next section](/docs/userguide/features/mfront/#preparing-opengeosys-to-use-new-mfront-model) for details on how to add new model).
The value of "parameter" attribute in `<material_property/>` has to match one of the content of tag `<name> </name>` in one of the [parameters](/docs/userguide/blocks/parameters/) block. For example:

```xml
<parameter>
    <name>ParameterName</name>
    <type>Constant</type>
    <value>value_of_parameter_Name</type>
</parameter>
```

Tag `<material_properties> </material_properties>` can be used as many times as it is necessary.
