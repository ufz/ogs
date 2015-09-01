#ifndef NUMLIB_INTERPOLATION_H
#define NUMLIB_INTERPOLATION_H

namespace NumLib
{

/**
 * Interpolates variables given at element nodes according to the given shape matrix.
 *
 * This function simply does the usual finite-element interpolation, i.e. multiplication
 * of nodal values with the shape function.
 *
 * @param nodal_values   vector of nodal values
 * @param shape_matrix_N shape matrix of the point to which will be interpolated
 * @param num_nodal_dof  number of nodal variables that will be interpolated
 * @param interpolated_values array of addresses to which the interpolated values will be written
 *
 * The size of interpolated_values must be equal to num_nodal_dof.
 */
template<typename NodalValues, typename ShapeMatrix>
void shapeFunctionInterpolate(
        const NodalValues& nodal_values,
        const ShapeMatrix& shape_matrix_N,
        const unsigned num_nodal_dof,
        double** interpolated_values
        )
{
    auto const num_nodes = shape_matrix_N.size();

    for (unsigned d=0; d<num_nodal_dof; ++d)
    {
        *interpolated_values[d] = 0.0;

        for (unsigned n=0; n<num_nodes; ++n)
        {
            *interpolated_values[d] += nodal_values[d*num_nodes+n] * shape_matrix_N(n);
        }
    }
}

}

#endif // NUMLIB_INTERPOLATION_H
