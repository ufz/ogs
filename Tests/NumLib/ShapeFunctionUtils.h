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

#include "NumLib/Fem/FiniteElement/ElementTraitsLagrange.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "Tests/Utils.h"

namespace ShapeFunctionUtils
{

template <unsigned GlobalDim, typename ElementType, typename PointsContainer>
auto computeShapeMatricesImpl(MeshLib::Element const& element,
                              PointsContainer const& points,
                              Type<ElementType>)
{
    using ET = NumLib::ElementTraitsLagrange<ElementType>;
    using ShpFct = typename ET::ShapeFunction;
    // Use dynamically allocated shape matrices s.t. their type does not depend
    // on the mesh element.
    using ShpMatPol = EigenDynamicShapeMatrixPolicy<ShpFct, GlobalDim>;

    bool const is_axially_symmetric = false;  // does not matter here
    return NumLib::computeShapeMatrices<ShpFct,
                                        ShpMatPol,
                                        GlobalDim,
                                        NumLib::ShapeMatrixType::N>(
        element, is_axially_symmetric, points);
}

template <typename BulkElementType, typename PointsContainer>
auto computeShapeMatricesBulk(BulkElementType const& e,
                              PointsContainer const& points)
{
    return computeShapeMatricesImpl<BulkElementType::dimension>(
        e, points, Type<BulkElementType>{});
}

template <typename PointsContainer>
auto computeShapeMatricesOnFace(MeshLib::Element const& face,
                                PointsContainer const& points)
{
    auto const compute =
        [&]<typename FaceElementType>(Type<FaceElementType> const type)
    {
        return computeShapeMatricesImpl<FaceElementType::dimension + 1>(
            face, points, type);
    };

    using namespace MeshLib;

    switch (auto const cell_type = face.getCellType())
    {
        case CellType::POINT1:
            return compute(Type<Point>{});
        case CellType::LINE2:
            return compute(Type<Line>{});
        case CellType::TRI3:
            return compute(Type<Tri>{});
        case CellType::QUAD4:
            return compute(Type<Quad>{});
        case CellType::TET4:
            return compute(Type<Tet>{});
        case CellType::HEX8:
            return compute(Type<Hex>{});
        case CellType::PRISM6:
            return compute(Type<Prism>{});
        case CellType::PYRAMID5:
            return compute(Type<Pyramid>{});
        default:
            OGS_FATAL("Unsupported cell type " + CellType2String(cell_type));
    }
}

// TODO return coord array
template <typename BulkElementType,
          typename ShpMat,
          unsigned GlobalDim = BulkElementType::dimension>
MathLib::Point3d interpolateNodeCoordinates(
    BulkElementType const& e,
    ShpMat const& N,
    [[maybe_unused]] std::integral_constant<unsigned, GlobalDim> global_dim =
        {})
{
    using ET = NumLib::ElementTraitsLagrange<BulkElementType>;
    using ShpFct = typename ET::ShapeFunction;
    using ShpMatPol = ShapeMatrixPolicyType<ShpFct, GlobalDim>;

    NumLib::TemplateIsoparametric<ShpFct, ShpMatPol> fe(e);
    return MathLib::Point3d(fe.interpolateCoordinates(N));
}

template <typename ShpMat>
auto interpolateNodeCoordinatesFace(MeshLib::Element const& face,
                                    ShpMat const& N)
{
    auto const interpolate =
        [&]<typename FaceElementType>(Type<FaceElementType>)
    {
        std::integral_constant<unsigned, FaceElementType::dimension + 1>
            global_dim{};
        return interpolateNodeCoordinates(
            static_cast<FaceElementType const&>(face), N, global_dim);
    };

    using namespace MeshLib;

    switch (auto const cell_type = face.getCellType())
    {
        case CellType::POINT1:
            return interpolate(Type<Point>{});
        case CellType::LINE2:
            return interpolate(Type<Line>{});
        case CellType::TRI3:
            return interpolate(Type<Tri>{});
        case CellType::QUAD4:
            return interpolate(Type<Quad>{});
        case CellType::TET4:
            return interpolate(Type<Tet>{});
        case CellType::HEX8:
            return interpolate(Type<Hex>{});
        case CellType::PRISM6:
            return interpolate(Type<Prism>{});
        case CellType::PYRAMID5:
            return interpolate(Type<Pyramid>{});
        default:
            OGS_FATAL("Unsupported cell type " + CellType2String(cell_type));
    }
}
}  // namespace ShapeFunctionUtils
