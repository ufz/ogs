/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_INTERPOLATION_H
#define NUMLIB_INTERPOLATION_H

#include<array>
#include<cassert>

namespace NumLib
{

/**
 * Interpolates variables given at element nodes according to the given shape matrix.
 *
 * This function simply does the usual finite-element interpolation, i.e. multiplication
 * of nodal values with the shape function.
 *
 * @param nodal_values   vector of nodal values, ordered by component
 * @param shape_matrix_N shape matrix of the point to which will be interpolated
 * @param interpolated_values array of addresses to which the interpolated values will be written
 */
template<typename NodalValues, typename ShapeMatrix, std::size_t NodalDOF>
void shapeFunctionInterpolate(
        const NodalValues& nodal_values,
        const ShapeMatrix& shape_matrix_N,
        std::array<double*, NodalDOF> interpolated_values
        )
{
    auto const num_nodes = shape_matrix_N.size();

    assert(num_nodes*NodalDOF == nodal_values.size());

    for (auto d=decltype(NodalDOF){0}; d<NodalDOF; ++d)
    {
        *interpolated_values[d] = 0.0;

        for (auto n=decltype(num_nodes){0}; n<num_nodes; ++n)
        {
            *interpolated_values[d] += nodal_values[d*num_nodes+n] * shape_matrix_N[n];
        }
    }
}

}

#endif // NUMLIB_INTERPOLATION_H
