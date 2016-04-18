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

namespace detail
{

template<unsigned DOFOffset, typename NodalValues, typename ShapeMatrix>
void shapeFunctionInterpolate(
        const NodalValues &/*nodal_values*/,
        const ShapeMatrix &/*shape_matrix_N*/)
{}

template<unsigned DOFOffset, typename NodalValues, typename ShapeMatrix, typename... Numbers>
void shapeFunctionInterpolate(
        const NodalValues &nodal_values,
        const ShapeMatrix &shape_matrix_N,
        double& interpolated_value,
        Numbers&... interpolated_values)
{
    auto const num_nodes = shape_matrix_N.size();
    double iv = 0.0;

    for (auto n=decltype(num_nodes){0}; n<num_nodes; ++n) {
        iv += nodal_values[DOFOffset*num_nodes+n] * shape_matrix_N[n];
    }

    interpolated_value = iv;

    shapeFunctionInterpolate<DOFOffset+1>(nodal_values, shape_matrix_N, interpolated_values...);
}

} // namespace detail

/*!
 * Interpolates variables given at element nodes according to the given shape matrix.
 *
 * This function simply does the usual finite-element interpolation, i.e. multiplication
 * of nodal values with the shape function.
 *
 * \param nodal_values vector of nodal values, ordered by component
 * \param shape_matrix_N shape matrix of the point to which will be interpolated
 * \param interpolated_value interpolated value of the first d.o.f. (output parameter)
 * \param interpolated_values interpolated value of further d.o.f. (output parameter)
 *
 * \note
 * \c nodal_values have to be ordered by component and it is assumed that all passed d.o.f. are
 * single-component and are interpolated using the same shape function.
 */
template<typename NodalValues, typename ShapeMatrix, typename... Numbers>
void shapeFunctionInterpolate(
        const NodalValues& nodal_values,
        const ShapeMatrix& shape_matrix_N,
        double& interpolated_value,
        Numbers&... interpolated_values
        )
{
    auto const num_nodal_dof = sizeof...(interpolated_values) + 1;
    auto const num_nodes = shape_matrix_N.size();

    assert(num_nodes*num_nodal_dof == nodal_values.size());
    (void) num_nodal_dof; // no warnings when not in debug build

    detail::shapeFunctionInterpolate<0>(nodal_values, shape_matrix_N, interpolated_value,
                                        interpolated_values...);
}

} // namespace NumLib

#endif // NUMLIB_INTERPOLATION_H
