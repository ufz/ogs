+++
date = "2018-02-27T11:00:13+01:00"
title = "Scaling source term"
author = "Feliks Kiszkurno"
weight = 13
+++

Consider a heat source defined as boundary in following mesh:

![Mesh](/docs/userguide/blocks/misc/figures/scalingsourceterm_mesh.png)
|:--:|
| *The "dent" in the right side of the mesh is a 2D, rectangular source* |

With height alongside Y-axis defined as $h=8$ meters and width alongside X-axis as $r=0.095$ meters.
Since the boundary parallel to Y-axis is [axially symmetrical](/docs/userguide/blocks/meshes/#axial-symmetry), the two dimensional rectangle will turn into a three dimensional cylinder after the mesh has been rotated by $360°$ with Y-axis boundary as rotation axis.

The source term has to take this into account.
Therefore its power has to be scaled - divided by the surface on which it will be delivered to the system.

Following this example. The surface of the cylinder is defined as follows:

$$A=2\pi r(h+r)$$

which will result in following surface:

$$A = 2\pi \cdot 0.095 \cdot (8 + 0.095) = 4.83$$

in $m^2$.

Let's assume that the heat supply is defined as 900 W.
Than we need to scale this value by the surface obtained from the previous equation:

$$S_{scaled}=900/4.83=186.34$$

and this is the value that has to be provided in the project file.

In this example the heat source is introduced as [Neumann boundary condition](/docs/userguide/blocks/boundary_conditions/#neumann) at part of the mesh. First, we can define a parameter called "heater" in the [parameters](/docs/userguide/blocks/parameters/) block:

```xml
<parameter>
    <name>heater</name>
    <!--values>900/(2*pi*0.095*(0.095+8)</values-->
    <!--values>186.34</values-->
    <type>Constant</type>
    <value>186.34</value>
</parameter>
```

Now we can refer to this value in the [process variables](/docs/userguide/blocks/process_variables/) block:

```xml
<process_variable>
    <name>temperature</name>
    <components>1</components>
    <order>1</order>
    <initial_condition>293.15</initial_condition>
    <boundary_conditions>
        <boundary_condition>
            <mesh>Mesh_Boundary</mesh>
            <type>Dirichlet</type>
            <component>0</component>
            <parameter>293.15</parameter>
        </boundary_condition>
        <boundary_condition>
            <mesh>Mesh_Heat_Source</mesh>
            <type>Neumann</type>
            <component>0</component>
            <parameter>heater</parameter>
        </boundary_condition>
    </boundary_conditions>
</process_variable>
```

Notice how different type of boundary conditions are used for the heat source and rest of the boundary.
More details on boundary conditions can be found [here](/docs/userguide/blocks/boundary_conditions/).
