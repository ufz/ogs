Type of mesh item for the data field to be set. Must be one of:

- `node` - for nodal data (point data in VTK)
- `cell` - for cell data (cell data in VTK)
- `integration_point` - for integration point data (field data in VTK). With a
  `<parameter_name>`, the parameter is evaluated once per element at the element
  centroid and applied to all of the element's integration points (see the
  `<set>` documentation).
