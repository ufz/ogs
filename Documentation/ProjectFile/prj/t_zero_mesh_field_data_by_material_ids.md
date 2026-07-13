\warning This tag has been removed. Use the `<overwrite_mesh_data>` element with a `<set>` (or `<take_original>`) sub-element and a `<material_ids>` entry instead.

It formerly specified material IDs for resetting the integration point data:
for the elements with the specified material IDs, their double-type integration
point data were set to zero. This is now expressed as a `<set>` operation (with
a zero-valued parameter) or handled explicitly per field.
