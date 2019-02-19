A local coordinate system used for space-dependent tensor transformations.

It offers a simple way for input of anisotropic tensors w.r.t. a coordinate
system.

The basis vectors form a transformation matrix \f$R = (e_0, e_1, e_2)\f$.

For a given anisotropic tensor \f$A\f$ parameter with the corresponding
\ref ogs_file_param__prj__parameters__parameter__use_local_coordinate_system
the tensor is rotated according to the formula: \f$A' = R\cdot A\cdot R^T\f$.

