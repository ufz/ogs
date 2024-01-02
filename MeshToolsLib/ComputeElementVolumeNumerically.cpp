/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 9, 2023, 2:01 PM
 */

#include "ComputeElementVolumeNumerically.h"

#include <typeinfo>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/MeshEnums.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/IntegrationMethodRegistry.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism15.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism6.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra13.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra5.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace MeshToolsLib
{
template <typename ShapeFunction>
double computeElementVolumeNumerically(MeshLib::Element const& e)
{
    // Space dimension is set to 3 in case that 1D or 2D element is inclined.
    constexpr int space_dim = 3;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, space_dim>;

    // Integration order is set to 3:
    auto const& integration_method =
        NumLib::IntegrationMethodRegistry::template getIntegrationMethod<
            typename ShapeFunction::MeshElement>(NumLib::IntegrationOrder{3});

    auto const shape_function_data =
        NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType, space_dim>(
            e, false /*is_axially_symmetric*/, integration_method);

    auto const n_integration_points = integration_method.getNumberOfPoints();
    double volume = 0.0;
    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const weight = integration_method.getWeightedPoint(ip).getWeight();
        volume += shape_function_data[ip].detJ * weight;
    }

    return volume;
}

double computeElementVolumeNumerically(MeshLib::Element const& e)
{
    switch (e.getCellType())
    {
        case MeshLib::CellType::LINE2:
            return computeElementVolumeNumerically<NumLib::ShapeLine2>(e);
        case MeshLib::CellType::LINE3:
            return computeElementVolumeNumerically<NumLib::ShapeLine3>(e);
        case MeshLib::CellType::TRI3:
            return computeElementVolumeNumerically<NumLib::ShapeTri3>(e);
        case MeshLib::CellType::TRI6:
            return computeElementVolumeNumerically<NumLib::ShapeTri6>(e);
        case MeshLib::CellType::QUAD4:
            return computeElementVolumeNumerically<NumLib::ShapeQuad4>(e);
        case MeshLib::CellType::QUAD8:
            return computeElementVolumeNumerically<NumLib::ShapeQuad8>(e);
        case MeshLib::CellType::QUAD9:
            return computeElementVolumeNumerically<NumLib::ShapeQuad9>(e);
        case MeshLib::CellType::TET4:
            return computeElementVolumeNumerically<NumLib::ShapeTet4>(e);
        case MeshLib::CellType::HEX8:
            return computeElementVolumeNumerically<NumLib::ShapeHex8>(e);
        case MeshLib::CellType::HEX20:
            return computeElementVolumeNumerically<NumLib::ShapeHex20>(e);
        case MeshLib::CellType::TET10:
            return computeElementVolumeNumerically<NumLib::ShapeTet10>(e);
        case MeshLib::CellType::PRISM6:
            return computeElementVolumeNumerically<NumLib::ShapePrism6>(e);
        case MeshLib::CellType::PRISM15:
            return computeElementVolumeNumerically<NumLib::ShapePrism15>(e);
        case MeshLib::CellType::PYRAMID5:
            return computeElementVolumeNumerically<NumLib::ShapePyra5>(e);
        case MeshLib::CellType::PYRAMID13:
            return computeElementVolumeNumerically<NumLib::ShapePyra13>(e);
        default:
            OGS_FATAL(
                "Numerical volume calculation is not available for element "
                "with type {}. ",
                MeshLib::CellType2String(e.getCellType()));
    }
}

}  // namespace MeshToolsLib
