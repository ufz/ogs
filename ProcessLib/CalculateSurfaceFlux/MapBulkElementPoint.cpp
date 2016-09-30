/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <array>

#include "MapBulkElementPoint.h"

namespace ProcessLib
{
MathLib::Point3d getBulkElementPoint(MeshLib::Quad const& quad,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint1D const& wp)
{
    switch (face_id) {
    case 0:
        return MathLib::Point3d{std::array<double, 3>{{wp[0], 0.0, 0.0}}};
    case 1:
        return MathLib::Point3d{std::array<double, 3>{{1.0, wp[0], 0.0}}};
    case 2:
        return MathLib::Point3d{std::array<double, 3>{{1.0 - wp[0], 1.0, 0.0}}};
    case 3:
        return MathLib::Point3d{std::array<double, 3>{{0.0, 1.0 - wp[0], 0.0}}};
    default:
        OGS_FATAL("Invalid face id '%u' for the quad.", face_id);
    }
}

MathLib::Point3d getBulkElementPoint(MeshLib::Hex const& hex,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint2D const& wp)
{
    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{
                std::array<double, 3>{{wp[1], wp[0], -1.0}}};
        case 1:
            return MathLib::Point3d{
                std::array<double, 3>{{wp[0], -1.0, wp[1]}}};
        case 2:
            return MathLib::Point3d{std::array<double, 3>{{1.0, wp[0], wp[1]}}};
        case 3:
            return MathLib::Point3d{
                std::array<double, 3>{{-wp[0], -1.0, wp[1]}}};
        case 4:
            return MathLib::Point3d{
                std::array<double, 3>{{-1.0, -wp[0], -wp[1]}}};
        case 5:
            return MathLib::Point3d{std::array<double, 3>{{wp[0], wp[1], 1.0}}};
        default:
            OGS_FATAL("Invalid face id '%u' for the hexahedron.", face_id);
    }
}

MathLib::Point3d getBulkElementPoint(MeshLib::Mesh const& mesh,
                                     std::size_t const bulk_element_id,
                                     std::size_t const bulk_face_id,
                                     MathLib::WeightedPoint1D const& wp)
{
    auto const* element = mesh.getElement(bulk_element_id);
    if (element->getCellType() == MeshLib::CellType::QUAD4)
    {
        MeshLib::Quad const& quad(*dynamic_cast<MeshLib::Quad const*>(element));
        return getBulkElementPoint(quad, bulk_face_id, wp);
    }
    OGS_FATAL("Wrong cell type '%s' or functionality not yet implemented.",
              MeshLib::CellType2String(element->getCellType()).c_str());
}

MathLib::Point3d getBulkElementPoint(MeshLib::Mesh const& mesh,
                                     std::size_t bulk_element_id,
                                     std::size_t bulk_face_id,
                                     MathLib::WeightedPoint2D const& wp)
{
    auto const* element = mesh.getElement(bulk_element_id);
    if (element->getCellType() == MeshLib::CellType::HEX8)
    {
        MeshLib::Hex const& hex(*dynamic_cast<MeshLib::Hex const*>(element));
        return getBulkElementPoint(hex, bulk_face_id, wp);
    }
    OGS_FATAL("Wrong cell type '%s' or functionality not yet implemented.",
              MeshLib::CellType2String(element->getCellType()).c_str());
}

// TODO disable the 3d elements in the local assembler creator
MathLib::Point3d getBulkElementPoint(MeshLib::Mesh const& mesh,
                                     std::size_t bulk_element_id,
                                     std::size_t bulk_face_id,
                                     MathLib::WeightedPoint3D const& wp)
{
    return MathLib::ORIGIN;
}
}  // end namespace ProcessLib
