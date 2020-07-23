/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include "MeshLib/Elements/Element.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"


namespace ProcessLib
{
template <typename ShapeFunction, typename ShapeMatricesType, int GlobalDim,
          NumLib::ShapeMatrixType SelectedShapeMatrixType =
              NumLib::ShapeMatrixType::ALL,
          typename IntegrationMethod>
std::vector<typename ShapeMatricesType::ShapeMatrices,
            Eigen::aligned_allocator<typename ShapeMatricesType::ShapeMatrices>>
initShapeMatrices(MeshLib::Element const& e, bool is_axially_symmetric,
                  IntegrationMethod const& integration_method)
{
    std::vector<
        typename ShapeMatricesType::ShapeMatrices,
        Eigen::aligned_allocator<typename ShapeMatricesType::ShapeMatrices>>
        shape_matrices;

    auto const fe =
        NumLib::createIsoparametricFiniteElement<ShapeFunction,
                                                 ShapeMatricesType>(e);

    unsigned const n_integration_points = integration_method.getNumberOfPoints();

    shape_matrices.reserve(n_integration_points);
    for (unsigned ip = 0; ip < n_integration_points; ++ip) {
        shape_matrices.emplace_back(ShapeFunction::DIM, GlobalDim,
                                     ShapeFunction::NPOINTS);
        fe.computeShapeFunctions(
                integration_method.getWeightedPoint(ip).getCoords(),
                shape_matrices[ip], GlobalDim, is_axially_symmetric);
    }

    return shape_matrices;
}

template <typename ShapeFunction, typename ShapeMatricesType>
double interpolateXCoordinate(
    MeshLib::Element const& e,
    typename ShapeMatricesType::ShapeMatrices::ShapeType const& N)
{
    auto const fe =
        NumLib::createIsoparametricFiniteElement<ShapeFunction,
                                                 ShapeMatricesType>(e);

    return fe.interpolateZerothCoordinate(N);
}

template <typename ShapeFunction, typename ShapeMatricesType>
std::array<double, 3> interpolateCoordinates(
    MeshLib::Element const& e,
    typename ShapeMatricesType::ShapeMatrices::ShapeType const& N)
{
    auto const fe =
        NumLib::createIsoparametricFiniteElement<ShapeFunction,
                                                 ShapeMatricesType>(e);

    return fe.interpolateCoordinates(N);
}

}  // namespace ProcessLib
