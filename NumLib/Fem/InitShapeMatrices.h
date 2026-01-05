// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <vector>

#include "FiniteElement/TemplateIsoparametric.h"
#include "MathLib/WeightedPoint.h"
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
            p.data(), shape_matrices.back(), GlobalDim, is_axially_symmetric);
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

    std::vector<MathLib::WeightedPoint> points;
    points.reserve(n_integration_points);
    for (int ip = 0; ip < n_integration_points; ++ip)
    {
        points.push_back(integration_method.getWeightedPoint(ip));
    }

    return computeShapeMatrices<ShapeFunction, ShapeMatricesType, GlobalDim,
                                SelectedShapeMatrixType>(
        e, is_axially_symmetric, points);
}

// Returned vector only contains one element.
template <typename ShapeFunction, typename ShapeMatricesType, int GlobalDim,
          ShapeMatrixType SelectedShapeMatrixType = ShapeMatrixType::ALL>
typename ShapeMatricesType::ShapeMatrices initShapeMatricesAtElementCenter(
    MeshLib::Element const& e, bool const is_axially_symmetric)
{
    static constexpr std::array<double, ShapeFunction::DIM> centre =
        ShapeFunction::reference_element_centre;

    static constexpr std::array integration_points = {
        MathLib::WeightedPoint{centre, 1.0}};

    auto const shape_matrices =
        computeShapeMatrices<ShapeFunction, ShapeMatricesType, GlobalDim,
                             SelectedShapeMatrixType>(e, is_axially_symmetric,
                                                      integration_points);
    return shape_matrices[0];
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
