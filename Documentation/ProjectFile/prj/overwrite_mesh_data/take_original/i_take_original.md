Restore the original mesh data for specified material IDs.

The `<take_original>` element is useful when nodes or elements along material boundaries have been overwritten by previous `<set>` operations and need to be restored to their original values.

**Required elements:**

- `<field_name>` - Name of the data field to restore
- `<mesh_item_type>` - Type of mesh item (node, cell, or `integration_point`)
- `<material_ids>` - Material IDs to restore (space-separated, supports ranges like `1:3` and `*`)

**Optional elements:**

- `<mesh>` - Name of the mesh (defaults to the first/bulk mesh)

**Example:**

```xml
<overwrite_mesh_data>
    <take_original>
        <field_name>temperature</field_name>
        <mesh_item_type>node</mesh_item_type>
        <material_ids>0</material_ids>
    </take_original>
</overwrite_mesh_data>
```
