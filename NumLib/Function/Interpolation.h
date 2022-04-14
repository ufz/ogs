/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>
#include <cassert>

#include "NumLib/Fem/CoordinatesMapping/NaturalNodeCoordinates.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace NumLib
{
namespace detail
{
//! \see ::NumLib::shapeFunctionInterpolate()
template <unsigned DOFOffset, typename NodalValues, typename ShapeMatrix>
void shapeFunctionInterpolate(const NodalValues& /*nodal_values*/,
                              const ShapeMatrix& /*shape_matrix_N*/)
{
}

//! \see ::NumLib::shapeFunctionInterpolate()
template <unsigned DOFOffset, typename NodalValues, typename ShapeMatrix,
          typename... ScalarTypes>
void shapeFunctionInterpolate(const NodalValues& nodal_values,
                              const ShapeMatrix& shape_matrix_N,
                              double& interpolated_value,
                              ScalarTypes&... interpolated_values)
{
    auto const num_nodes = shape_matrix_N.size();
    double iv = 0.0;

    for (auto n = decltype(num_nodes){0}; n < num_nodes; ++n)
    {
        iv += nodal_values[DOFOffset * num_nodes + n] * shape_matrix_N[n];
    }

    interpolated_value = iv;

    shapeFunctionInterpolate<DOFOffset + 1>(nodal_values, shape_matrix_N,
                                            interpolated_values...);
}

}  // namespace detail

/*!
 * Interpolates variables given at element nodes according to the given shape
 * matrix.
 *
 * This function simply does the usual finite-element interpolation, i.e.
 * multiplication of nodal values with the shape function.
 *
 * \param nodal_values vector of nodal values, ordered by component
 * \param shape_matrix_N shape matrix of the point to which will be interpolated
 * \param interpolated_value interpolated value of the first d.o.f. (output
 * parameter) \param interpolated_values interpolated value of further d.o.f.
 * (output parameter)
 *
 * \tparam NodalValues  type of the container where nodal values are stored
 * \tparam ShapeMatrix  type of the shape matrix \f$N\f$.
 * \tparam ScalarValues all of the types in this pack must currently be \c
 * double.
 *
 * \note
 * \c nodal_values have to be ordered by component and it is assumed that all
 * passed d.o.f. are single-component and are interpolated using the same shape
 * function.
 */
template <typename NodalValues, typename ShapeMatrix, typename... ScalarTypes>
void shapeFunctionInterpolate(const NodalValues& nodal_values,
                              const ShapeMatrix& shape_matrix_N,
                              double& interpolated_value,
                              ScalarTypes&... interpolated_values)
{
    auto const num_nodal_dof = sizeof...(interpolated_values) + 1;
    auto const num_nodes = shape_matrix_N.size();

    assert(num_nodes * num_nodal_dof ==
           static_cast<std::size_t>(nodal_values.size()));
    (void)num_nodal_dof;
    (void)num_nodes;  // no warnings when not in debug build

    detail::shapeFunctionInterpolate<0>(nodal_values, shape_matrix_N,
                                        interpolated_value,
                                        interpolated_values...);
}

/// Interpolates scalar \c node_values given in lower order element nodes (e.g.
/// the base nodes) to higher order element's nodes (e.g. quadratic nodes) and
/// writes the result into the global property vector.
///
/// The base nodes' values are copied.  For each higher order node the shape
/// matrices are evaluated for the lower order element (the base nodes), and
/// used for the the scalar quantity interpolation.
template <typename LowerOrderShapeFunction, typename HigherOrderMeshElementType,
          int GlobalDim, typename EigenMatrixType>
void interpolateToHigherOrderNodes(
    MeshLib::Element const& element, bool const is_axially_symmetric,
    Eigen::MatrixBase<EigenMatrixType> const& node_values,
    MeshLib::PropertyVector<double>& interpolated_values_global_vector)
{
    assert(dynamic_cast<HigherOrderMeshElementType const*>(&element));
    assert(node_values.cols() == 1);  // Scalar quantity only.

    using SF = LowerOrderShapeFunction;
    using ShapeMatricesType = ShapeMatrixPolicyType<SF, GlobalDim>;

    int const number_base_nodes = element.getNumberOfBaseNodes();
    int const number_all_nodes = element.getNumberOfNodes();

    // Copy the values for linear nodes.
    for (int n = 0; n < number_base_nodes; ++n)
    {
        std::size_t const global_index = getNodeIndex(element, n);
        interpolated_values_global_vector[global_index] = node_values[n];
    }

    //
    // Interpolate values for higher order nodes.
    //

    int const number_higher_order_nodes = number_all_nodes - number_base_nodes;
    std::vector<MathLib::Point3d> higher_order_nodes;
    higher_order_nodes.reserve(number_higher_order_nodes);
    for (int n = 0; n < number_higher_order_nodes; ++n)
    {
        higher_order_nodes.emplace_back(
            NaturalCoordinates<HigherOrderMeshElementType>::coordinates
                [number_base_nodes + n]);
    }

    // Shape matrices evaluated at higher order nodes' coordinates.
    auto const shape_matrices =
        computeShapeMatrices<SF, ShapeMatricesType, GlobalDim,
                             ShapeMatrixType::N>(element, is_axially_symmetric,
                                                 higher_order_nodes);

    for (int n = 0; n < number_higher_order_nodes; ++n)
    {
        std::size_t const global_index =
            getNodeIndex(element, number_base_nodes + n);
        interpolated_values_global_vector[global_index] =
            shape_matrices[n].N * node_values;
    }
}

}  // namespace NumLib
