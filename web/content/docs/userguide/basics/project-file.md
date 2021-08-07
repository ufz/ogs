+++
date = "2018-11-14T11:00:13+01`:00"
title = "Project file syntax"
author = "Lars Bilke"
weight = 4

[menu]
  [menu.userguide]
    parent = "basics"
+++

Project files `.prj` have to be valid XML documents. For information of specific tags see our [Doxygen-documentation](https://doxygen.opengeosys.org/d1/d91/ogs_file_param__projectfile).

Two methods allow you to modularize your project files and avoid repetition:

## Option 1: Include XML-content from other files

The `<include file="other_file.xml" />`-tag allows to include the content of another XML into the current project file.

**Limitation:** Only one child `include`-element per (regular) XML-element is allowed.

### Example

`Tests/Data/Elliptic/circle_radius_1/circle_1e1_axi.prj`:

```xml
    <processes>
        <include file="../cube_1x1x1_SteadyStateDiffusion/SteadyStateDiffusion.xml"/>
    </processes>
```

`Tests/Data/Elliptic/cube_1x1x1_SteadyStateDiffusion/SteadyStateDiffusion.xml`:
```xml
<process>
    <name>SteadyStateDiffusion</name>
    <type>STEADY_STATE_DIFFUSION</type>
    <integration_order>2</integration_order>
    <process_variables>
        <process_variable>pressure</process_variable>
    </process_variables>
    <secondary_variables>
        <secondary_variable internal_name="darcy_velocity" output_name="v"/>
    </secondary_variables>
</process>
```

## Option 2: Apply patch files to the project file

Patch files contain `<replace>`, `<add>` and `<remove>`-elements with [XPath](https://en.wikipedia.org/wiki/XPath)-selectors to modify a specific part of the project file (in-memory during run-time):

`Tests/Data/Elliptic/square_1x1_SteadyStateDiffusion/square_neumann.xml`:
```xml
<?xml version="1.0" encoding="ISO-8859-1"?>
<OpenGeoSysProjectDiff>
    <add sel="/*/time_loop/output/prefix/text()" pos="after">_neumann</add>
    <replace sel="/*/parameters/parameter[2]/name/text()">p_neumann</replace>
    <replace sel="/*/parameters/parameter[3]/name/text()">p_Dirichlet</replace>
    <replace sel="/*/parameters/parameter[3]/value/text()">1</replace>
    <replace sel="/*/process_variables/process_variable/boundary_conditions/boundary_condition[1]/parameter/text()">p_Dirichlet</replace>
    <replace sel="/*/process_variables/process_variable/boundary_conditions/boundary_condition[2]/type/text()">Neumann</replace>
    <replace sel="/*/process_variables/process_variable/boundary_conditions/boundary_condition[2]/parameter/text()">p_neumann</replace>
    <add sel="/*/process_variables/process_variable/boundary_conditions/boundary_condition[1]" pos="after">
        <boundary_condition>
            <geometrical_set>square_1x1_geometry</geometrical_set>
            <geometry>bottom</geometry>
            <type>Dirichlet</type>
            <parameter>p_Dirichlet</parameter>
        </boundary_condition>
    </add>
</OpenGeoSysProjectDiff>
```

### Option 2a: Supply patch file(s) additionally to the prj-file

```bash
ogs -p path/to/square_1e0_neumann.xml [other/optional/patch_file.xml] path/to/square_1e0.prj
```

Supplied patch files are applied in the given order.

### Option 2b: Use a patch file directly

If the patch file specifies a `base_file`:

`Tests/Data/Elliptic/cube_1x1x1_SteadyStateDiffusion/cube_1e1.xml`:
```xml
<?xml version="1.0" encoding="ISO-8859-1"?>
<OpenGeoSysProjectDiff base_file="cube_1e0.prj">
    <replace sel="/*/mesh/text()">cube_1x1x1_hex_1e1.vtu</replace>
    <replace sel="/*/time_loop/output/prefix/text()">cube_1e1</replace>
</OpenGeoSysProjectDiff>

```

you can just pass this patch file:

```bash
ogs path/to/cube_1e1.xml
```

In this case you have just one patch file.

---

Both methods, include and patch, can be combined.

---

## Check project file syntax with `xmllint`

 You can check the formatting with the [`xmllint`-tool](http://xmlsoft.org/xmllint.html):

```bash
xmllint --noout myproj.prj
```

### Install `xmllint`

<div class='win'>

We recommend to install via [Chocolatey](https://chocolatey.org):

```powershell
choco install xsltproc
```

<div class='note'>

### <i class="far fa-info-circle"></i> Alternative installation

Another method is to use the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) where you can simply install Linux packages:

```bash
sudo apt-get install libxml2-utils
```

</div>

</div>

<div class='linux'>

```bash
sudo apt-get install libxml2-utils
```

</div>

<div class='mac'>
`xmllint` is part of the Homebrew xmlstarlet package:

```bash
brew install xmlstarlet
```

</div>
