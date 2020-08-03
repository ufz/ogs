+++
date = "2020-02-11T10:44:18+01:00"
title = "Linear; Element deactivation"
project = "Mechanics/Linear/square_with_deactivated_hole.prj"
author = "Wenqing Wang"
weight = 111

[menu]
  [menu.benchmarks]
    parent = "small-deformations"

+++

{{< data-link >}}

The problem definition is exactly the same as that of [Disc with hole]({{<ref "mechanics-linear-disc-with-hole.md">}}). With the element deactivation approach, the problem is solved as 2D and 3D benchmarks, respectively.

The input data set of the element deactivation approach is specified inside the tag of   `<process_variable> ... </process_variable>`. For example, the following input means to deactivate the elements with MaterialIDs of 0  within a time interval from 0 to 1:

```xml
    <process_variables>
        <process_variable>
            ...
            <deactivated_subdomains>
                <deactivated_subdomain>
                    <time_interval>
                        <start>0</start>
                        <end>1</end>
                    </time_interval>
                    <material_ids>0</material_ids>
                </deactivated_subdomain>
            </deactivated_subdomains>
            ....
        <process_variables>
    <process_variable>
```

 The input syntax can be also found in the project files of the two benchmarks:

* [2D](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/Tests/Data/Mechanics/Linear/square_with_deactivated_hole.prj)
* [3D](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/Tests/Data/Mechanics/Linear/ElementDeactivation3D/element_deactivation_M_3D.prj)

## Mesh

![2D and 3D meshes](../element_deactivation_2D_3D_mesh.png)

## Results and evaluation

![2D results](../element_deactivation_2D.png)

![3D results:](../element_deactivation_3D.png)
