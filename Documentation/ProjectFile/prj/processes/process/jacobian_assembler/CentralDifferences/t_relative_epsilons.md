Specifies the magnitudes of the perturbations used to compute the numerical
Jacobian.

**Note**: The deformation block of the local Jacobian matrix has already been
 calculated analytically. Therefore, the perturbations for the deformation
 components is not needed for computation. However, they are still required as
 input and can be set to zeros.

The magnitudes are specified relative to the \c `component_magnitudes`.
The number of values given must match the one of the \c `component_magnitudes`.
