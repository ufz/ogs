/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include "FiniteElement/TemplateIsoparametric.h"
#include "MeshLib/Elements/Element.h"

namespace NumLib
{
template <typename ShapeFunction, typename ShapeMatricesType, int GlobalDim,
          ShapeMatrixType SelectedShapeMatrixType = ShapeMatrixType::ALL,
          typename PointContainer>
std::vector<typename ShapeMatricesType::ShapeMatrices,
            Eigen::aligned_allocator<typename ShapeMatricesType::ShapeMatrices>>
computeShapeMatrices(MeshLib::Element const& e, bool const is_axially_symmetric,
                     PointContainer const& points)
{
    std::vector<
        typename ShapeMatricesType::ShapeMatrices,
        Eigen::aligned_allocator<typename ShapeMatricesType::ShapeMatrices>>
        shape_matrices;

    auto const fe =
        createIsoparametricFiniteElement<ShapeFunction, ShapeMatricesType>(e);

    shape_matrices.reserve(points.size());
    for (auto const& p : points)
    {
        shape_matrices.emplace_back(ShapeFunction::DIM, GlobalDim,
                                    ShapeFunction::NPOINTS);
        fe.template computeShapeFunctions<SelectedShapeMatrixType>(
            p.getCoords(), shape_matrices.back(), GlobalDim,
            is_axially_symmetric);
    }

    return shape_matrices;
}

template <typename ShapeFunction, typename ShapeMatricesType, int GlobalDim,
          ShapeMatrixType SelectedShapeMatrixType = ShapeMatrixType::ALL,
          typename IntegrationMethod>
std::vector<typename ShapeMatricesType::ShapeMatrices,
            Eigen::aligned_allocator<typename ShapeMatricesType::ShapeMatrices>>
initShapeMatrices(MeshLib::Element const& e, bool const is_axially_symmetric,
                  IntegrationMethod const& integration_method)
{
    int const n_integration_points = integration_method.getNumberOfPoints();

    std::vector<typename IntegrationMethod::WeightedPoint> points;
    points.reserve(n_integration_points);
    for (int ip = 0; ip < n_integration_points; ++ip)
    {
        points.push_back(integration_method.getWeightedPoint(ip));
    }

    return computeShapeMatrices<ShapeFunction, ShapeMatricesType, GlobalDim,
                                SelectedShapeMatrixType>(
        e, is_axially_symmetric, points);
}

template <typename ShapeFunction, typename ShapeMatricesType>
double interpolateXCoordinate(
    MeshLib::Element const& e,
    typename ShapeMatricesType::ShapeMatrices::ShapeType const& N)
{
    auto const fe =
        createIsoparametricFiniteElement<ShapeFunction, ShapeMatricesType>(e);

    return fe.interpolateZerothCoordinate(N);
}

template <typename ShapeFunction, typename ShapeMatricesType>
std::array<double, 3> interpolateCoordinates(
    MeshLib::Element const& e,
    typename ShapeMatricesType::ShapeMatrices::ShapeType const& N)
{
    auto const fe =
        createIsoparametricFiniteElement<ShapeFunction, ShapeMatricesType>(e);

    return fe.interpolateCoordinates(N);
}

}  // namespace NumLib
