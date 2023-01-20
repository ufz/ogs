This is a link to the name of the geometrical object within the geometrical
set specified in the `<name>` tag.

In the example gml snippet

```xml
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>cube_geometry</name>
    <points>
    <!-- definition of points -->
    </points>

    <surfaces>
        <surface id="0" name="leftsurface">
            <element p1="0" p2="1" p3="2"/>
            <element p1="0" p2="3" p3="2"/>
        </surface>
    </surfaces>
</OpenGeoSysGLI>
```

the name of the geometry is __leftsurface__.
