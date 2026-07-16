Set a mesh data field from a parameter value for specified material IDs.

The `<set>` element assigns values from a parameter to the specified mesh data field for elements/nodes belonging to the given material IDs.

**Required elements:**

- `<field_name>` - Name of the data field to set
- `<mesh_item_type>` - Type of mesh item (node, cell, or `integration_point`)
- `<material_ids>` - Material IDs to apply the data to (space-separated, supports ranges like `1:3` and `*`)
- `<parameter_name>` - Name of the parameter providing the values

**Optional elements:**

- `<mesh>` - Name of the mesh (defaults to the first/bulk mesh)

**Note on `integration_point` data:** the parameter is evaluated once per
element, at the element centroid, and that single value is written to every
integration point of the element. This is exact for the intended use case of a
constant value across a material; a spatially varying parameter is not resolved
per integration point.

**Example:**

```xml
<overwrite_mesh_data>
    <set>
        <field_name>temperature</field_name>
        <mesh_item_type>node</mesh_item_type>
        <material_ids>1</material_ids>
        <parameter_name>temp_param</parameter_name>
    </set>
</overwrite_mesh_data>
```
