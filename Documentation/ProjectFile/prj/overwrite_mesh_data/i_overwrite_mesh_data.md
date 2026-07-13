Overwrite mesh data based on material ids.

The `<overwrite_mesh_data>` tag contains `<remove>`, `<set>`, and/or `<take_original>` elements that define how to modify mesh data fields for specific material IDs or mesh items.

**Available operations:**

- `<remove>` - Remove a mesh data field entirely
- `<set>` - Set a mesh data field from a parameter value
- `<take_original>` - Restore the original mesh data (useful after previous overwrites)

**Example:**

```xml
<overwrite_mesh_data>
    <remove>
        <field_name>old_data</field_name>
        <mesh_item_type>cell</mesh_item_type>
    </remove>
    <set>
        <field_name>temperature</field_name>
        <mesh_item_type>node</mesh_item_type>
        <material_ids>1</material_ids>
        <parameter_name>temp_param</parameter_name>
    </set>
    <take_original>
        <field_name>temperature</field_name>
        <mesh_item_type>node</mesh_item_type>
        <material_ids>0</material_ids>
    </take_original>
</overwrite_mesh_data>
```
