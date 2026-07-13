Remove the specified mesh data field from the mesh.

The `<field_name>` specifies which field to remove, and `<mesh_item_type>` specifies whether it is a node, cell, or integration point data field.

**Required elements:**

- `<field_name>` - Name of the data field to remove
- `<mesh_item_type>` - Type of mesh item (node, cell, or `integration_point`)

**Optional elements:**

- `<mesh>` - Name of the mesh (defaults to the first/bulk mesh)

**Example:**

```xml
<overwrite_mesh_data>
    <remove>
        <field_name>old_data</field_name>
        <mesh_item_type>cell</mesh_item_type>
    </remove>
</overwrite_mesh_data>
```
