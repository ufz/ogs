+++
date = "2023-03-09 09:22:24"
title = "Thermal conductivity: effective porosity mixing"
author = "Feliks Kiszkurno"
weight = 6
+++

In some process, the effective thermal conductivity can be calculated automatically depending on the conductivities of solid and liquid phases and porosity.

Following example can be considered: the Layer0 is a porous clay fully saturated by water. In such a case, in order to run simulation for correct value of thermal conductivity for such a medium, there has to be separate values for thermal conductivity of water (in liquid phase) and clay (in solid phase) defined. Those two values together with porosity can be used to obtain parameter representative for the whole medium. It can be done with following equation for volumetric mixing:

$$
\lambda_{medium}=\lambda_{water}*\phi+\lambda_{clay}\cdot(1-\phi)
$$

where $\lambda$ indicates thermal conductivity and $\phi$ indicates porosity. OpenGeoSys can do this internally. The requirement for it to work is that both phases have property with ```<name>thermal_conductivity</name>``` and porosity is defined for the whole medium. Than $\lambda$ for the whole medium can be defined as follows:

```xml
<property>
    <name>thermal_conductivity</name>
    <type>EffectiveThermalConductivityPorosityMixing</type>
</property>
```

This is how a media block with all required elements to use thermal porosity mixing can be defined:

```xml
<medium>
    <phases>
        <phase>
            <type>AqueousLiquid</type>
            <properties>
                <property>
                    <name>thermal_conductivity</name>
                    <type>Constant</type>
                    <value>thermal_conductivity_liquid_value</value>
                </property>
            </properties>
        <phase>
        <phase>
            <type>Solid</type>
            <properties>
                <property>
                    <name>thermal_conductivity</name>
                    <type>Constant</type>
                    <value>thermal_conductivity_solid_value</value>
                </property>
            </properties>
        </phase>
    </phases>
    <properties>
        <property>
            <name>porosity</name>
            <type>Constant</type>
            <value>medium_porosity_value</value>
        </property>
        <property>
            <name>thermal_conductivity</name>
            <type>EffectiveThermalConductivityPorosityMixing</type>
        </property>
    </properties>
</medium>
```
