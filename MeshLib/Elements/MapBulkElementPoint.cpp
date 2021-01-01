/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <array>

#include "MapBulkElementPoint.h"

namespace MeshLib
{
MathLib::Point3d getBulkElementPoint(MeshLib::Tri const& /*tri*/,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint1D const& wp)
{
    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{std::array<double, 3>{{wp[0], 0.0, 0.0}}};
        case 1:
            return MathLib::Point3d{
                std::array<double, 3>{{1 - wp[0], wp[0], 0.0}}};
        case 2:
            return MathLib::Point3d{
                std::array<double, 3>{{0.0, 1 - wp[0], 0.0}}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the tri.", face_id);
    }
}

MathLib::Point3d getBulkElementPoint(MeshLib::Quad const& /*quad*/,
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
        OGS_FATAL("Invalid face id '{:d}' for the quad.", face_id);
    }
}

MathLib::Point3d getBulkElementPoint(MeshLib::Hex const& /*hex*/,
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
            OGS_FATAL("Invalid face id '{:d}' for the hexahedron.", face_id);
    }
}

MathLib::Point3d getBulkElementPoint(MeshLib::Prism const& /*prism*/,
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
                std::array<double, 3>{{wp[0]/2.0 + 0.5, 0.0, wp[1]}}};
        case 2:
            return MathLib::Point3d{
                std::array<double, 3>{{0.5 - wp[0]/2.0, 0.5 + wp[0]/2.0, wp[1]}}};
        case 3:
            return MathLib::Point3d{
                std::array<double, 3>{{0, -wp[0]/2.0 + 0.5, wp[1]}}};
        case 4:
            return MathLib::Point3d{
                std::array<double, 3>{{wp[0], wp[1], 1.0}}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the prism.", face_id);
    }
}

MathLib::Point3d getBulkElementPoint(MeshLib::Pyramid const& /*pyramid*/,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint2D const& wp)
{
    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{
                std::array<double, 3>{{2 * wp[0] - 1, -1.0, 2 * wp[1] - 1}}};
        case 1:
            return MathLib::Point3d{
                std::array<double, 3>{{1.0, 2 * wp[0] - 1, 2 * wp[1] - 1}}};
        case 2:
            return MathLib::Point3d{
                std::array<double, 3>{{1 - 2 * wp[0], 1.0, 2 * wp[1] - 1}}};
        case 3:
            return MathLib::Point3d{std::array<double, 3>{
                {-1.0, 2 * wp[1] - 1, 2 * wp[0] - 1}}};
        case 4:
            return MathLib::Point3d{
                std::array<double, 3>{{-wp[0], wp[1], -1.0}}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the pyramid.", face_id);
    }
}

MathLib::Point3d getBulkElementPoint(MeshLib::Tet const& /*tet*/,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint2D const& wp)
{
    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{std::array<double, 3>{{wp[1], wp[0], 0.0}}};
        case 1:
            return MathLib::Point3d{std::array<double, 3>{{wp[0], 0.0, wp[1]}}};
        case 2:
            return MathLib::Point3d{
                std::array<double, 3>{{1 - wp[0] - wp[1], wp[0], wp[1]}}};
        case 3:
            return MathLib::Point3d{
                std::array<double, 3>{{0, wp[1], wp[0]}}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the tetrahedron.", face_id);
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
    if (element->getCellType() == MeshLib::CellType::TRI3)
    {
        MeshLib::Tri const& tri = *static_cast<MeshLib::Tri const*>(element);
        return getBulkElementPoint(tri, bulk_face_id, wp);
    }
    OGS_FATAL("Wrong cell type '{:s}' or functionality not yet implemented.",
              MeshLib::CellType2String(element->getCellType()));
}

MathLib::Point3d getBulkElementPoint(MeshLib::Mesh const& mesh,
                                     std::size_t bulk_element_id,
                                     std::size_t bulk_face_id,
                                     MathLib::WeightedPoint2D const& wp)
{
    auto const* element = mesh.getElement(bulk_element_id);
    if (element->getCellType() == MeshLib::CellType::HEX8)
    {
        MeshLib::Hex const& hex = *static_cast<MeshLib::Hex const*>(element);
        return getBulkElementPoint(hex, bulk_face_id, wp);
    }
    if (element->getCellType() == MeshLib::CellType::PRISM6)
    {
        MeshLib::Prism const& prism =
            *static_cast<MeshLib::Prism const*>(element);
        return getBulkElementPoint(prism, bulk_face_id, wp);
    }
    if (element->getCellType() == MeshLib::CellType::PYRAMID5)
    {
        MeshLib::Pyramid const& pyramid =
            *static_cast<MeshLib::Pyramid const*>(element);
        return getBulkElementPoint(pyramid, bulk_face_id, wp);
    }
    if (element->getCellType() == MeshLib::CellType::TET4)
    {
        MeshLib::Tet const& tet = *static_cast<MeshLib::Tet const*>(element);
        return getBulkElementPoint(tet, bulk_face_id, wp);
    }
    OGS_FATAL("Wrong cell type '{:s}' or functionality not yet implemented.",
              MeshLib::CellType2String(element->getCellType()));
}

// TODO disable the 3d elements in the local assembler creator
MathLib::Point3d getBulkElementPoint(MeshLib::Mesh const& /*mesh*/,
                                     std::size_t /*bulk_element_id*/,
                                     std::size_t /*bulk_face_id*/,
                                     MathLib::WeightedPoint3D const& /*wp*/)
{
    return MathLib::ORIGIN;
}
}  // namespace MeshLib
