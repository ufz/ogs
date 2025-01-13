+++
date = "2018-02-27T11:00:13+01:00"
title = "Project file structure"
author = "Feliks Kiszkurno"
weight = 31
+++

The main input taken by OpenGeoSys is a path to a project file. This file has a ".prj" extension and describes all aspects of
the simulation that OpenGeoSys is expected to execute.

The syntax of the project file follows the syntax of a typical XML file and can be edited with any editor capable of working
with this kind of files.

The structure of the project file is relatively flexible: as long as the semantics of the file is intact, the main blocks can be
provided in any order.

Depending on the scope of the simulation, the content of the project file will differ. However, there are some basic requirement that have to be fulfilled.

<div class="note">
Please note, that this page is not a complete description of the XML standard. It only provides some information to interact
with and understand OpenGeoSys "*.prj" files. More meta-information on XML can be found here: <a href="https://www.w3.org/standards/xml/">XML Standard</a>. <!-- TODO: Consider giving a direct link to an xml-tutorial. -->
</div>

The XML files structure follows an abstract idea of a tree: i.e., here one starts from the roots and continues up to the
leaves. As such the main tag within an XML file is called root. In OpenGeoSys project files the root tag looks like this: `<
OpenGeoSysProject> </OpenGeoSysProject>`. Between opening and closing tags all the project relevant information is contained.
Everything between those tags follows a hierarchical structure. There will be thematic blocks that cover specific aspects of
the simulation. They define for example mesh files, or a process used in the simulation. The following XML-snippet shows a
small structure containing the two elementary blocks, which are empty for now:

```xml
<OpenGeoSysProject>
    <mesh>
    ...
    </mesh>

    <process>
    ...
    </process>
</OpenGeoSysProject>
```

The blocks `<mesh> </mesh>` and `<processes> </processes>` are at the same "level", which is indicated by the indent. Those are
the child elements of the root `<OpenGeoSysProject>`, which is their parent element. As they are on the same level, they can be
described as siblings. This hierarchy is recursive. This means, that each of the child elements can have it's own child
elements and be parent to them (think of family structure: grandparents -> parent -> children, where the children are
grandchildren of the grandparents).

In following snippet you can see a `<processes> </processes>` element containing child element `<processes> </processes>`, which again contains two child elements: `<name> </name>` and `<type> </type>`:

```xml
<OpenGeoSysProject>
    <process>
        <process>
            <name>SteadyStateDiffusion</name>
            <type>STEADY_STATE_DIFFUSION</type>
        </process>
    </process>
</OpenGeoSysProject>
```

This snippet can be understood in more human readable way as follows:

*In this project file used by OpeGeoSys there is one instance of `<processes>`, which we call SteadyStateDiffusion and this
 instance is a realization of the type STEADY_STATE_DIFFUSION implemented in OpenGeoSys.*

## Attributes

If a value is set between `<` and `>`, beside the tag, it is called attribute. For example:

```xml
<planet satellite_name="Moon">
```

The entity `planet` has the attribute `satellite_name`, which is set to `Moon`.

In project files, attributes are used in some blocks to distinguished between multiple instances of the same entity. For
example [media](/docs/userguide/blocks/media/#media) are identified by their attribute "id":

```xml
<medium id="0">
```

and processes in [time loop](/docs/userguide/blocks/time_loop/#process) by the "ref" attribute:

```xml
<process ref="SteadyStateDiffusion">
```

## Comments

If you want to switch between two version of a certain part of the project file, two values etc, but without need to remove those fragments completely from the file, you can comment them out:

```xml
<!--value>this value will be comment out and ignored by OpenGeoSys</value-->
<value>this value will be read and used by OpenGeoSys</value>
```

Please note that comments cannot be nested. If you try to comment out a block that already contains a commented section, it will result in an error.

## More information

More details on what blocks are required/available in the project file and what they can contain, please see: [Project file over view](/docs/userguide/blocks/intro/).

## File encoding

In order to ensure correct encoding, it is good practice to add following line at the beginning of the project file:

```xml
<?xml version="1.0" encoding="ISO-8859-1"?>
```
