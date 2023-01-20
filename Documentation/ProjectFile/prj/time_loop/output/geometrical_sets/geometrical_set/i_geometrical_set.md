Tag that describes the mesh output according to a geometric object (point, polyline,
(plane) surface). The geometric object is defined in a gml file.

```xml
<geometrical_sets>
    <geometrical_set>
        <name>cube_geometry</name>
        <geometry>leftsurface</geometry>
    </geometrical_set>
    ...
    <geometrical_set>
        <name>cube_geometry</name>
        <geometry>rightsurface</geometry>
    </geometrical_set>
</geometrical_sets>
```

The name of the internal generated mesh is constructed from the name of the
geometrical set, an underscore, and the name of the geometrical object.
In the case above the internal mesh name would be __cube_geometry_leftsurface__.
