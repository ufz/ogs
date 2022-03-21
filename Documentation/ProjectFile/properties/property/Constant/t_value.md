The numerical value of one of the implemented data types, as given in MaterialPropertyLib::PropertyDataType.

Scalars, vectors (2 or 3 components), symmetric 3D tensors (6 components) and
full tensors (4 or 9 components) are supported. Based on the number of
components, OGS will deduce the type of the quantity (scalar/vector/tensor,
1D/2D/3D) automatically.

The elements of `<values>` are in row-major order. I.e., if you write the
following into your project file
```
<value>
1 2
3 4
</value>
```
OGS will read the proper 2x2 matrix.

For the order of symmetric tensor components please refer to [the user guide](https://www.opengeosys.org/docs/userguide/basics/conventions/#symmetric-tensors).
