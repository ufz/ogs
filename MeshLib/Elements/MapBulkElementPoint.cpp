/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MapBulkElementPoint.h"

#include <array>

namespace MeshLib
{
MathLib::Point3d getBulkElementPoint(Line const& /*line*/,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& /*wp*/)
{
    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{std::array<double, 3>{{-1.0, 0.0, 0.0}}};
        case 1:
            return MathLib::Point3d{std::array<double, 3>{{1.0, 0.0, 0.0}}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for a line.", face_id);
    }
}

MathLib::Point3d getBulkElementPoint(Tri const& /*tri*/,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp)
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

MathLib::Point3d getBulkElementPoint(Quad const& /*quad*/,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp)
{
    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{std::array<double, 3>{{wp[0], 0.0, 0.0}}};
        case 1:
            return MathLib::Point3d{std::array<double, 3>{{1.0, wp[0], 0.0}}};
        case 2:
            return MathLib::Point3d{
                std::array<double, 3>{{1.0 - wp[0], 1.0, 0.0}}};
        case 3:
            return MathLib::Point3d{
                std::array<double, 3>{{0.0, 1.0 - wp[0], 0.0}}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the quad.", face_id);
    }
}

MathLib::Point3d getBulkElementPoint(Hex const& /*hex*/,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp)
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

MathLib::Point3d getBulkElementPoint(Prism const& /*prism*/,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp)
{
    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{
                std::array<double, 3>{{wp[1], wp[0], -1.0}}};
        case 1:
            return MathLib::Point3d{
                std::array<double, 3>{{wp[0] / 2.0 + 0.5, 0.0, wp[1]}}};
        case 2:
            return MathLib::Point3d{std::array<double, 3>{
                {0.5 - wp[0] / 2.0, 0.5 + wp[0] / 2.0, wp[1]}}};
        case 3:
            return MathLib::Point3d{
                std::array<double, 3>{{0, -wp[0] / 2.0 + 0.5, wp[1]}}};
        case 4:
            return MathLib::Point3d{std::array<double, 3>{{wp[0], wp[1], 1.0}}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the prism.", face_id);
    }
}

MathLib::Point3d getBulkElementPoint(Pyramid const& /*pyramid*/,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp)
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
            return MathLib::Point3d{
                std::array<double, 3>{{-1.0, 2 * wp[1] - 1, 2 * wp[0] - 1}}};
        case 4:
            return MathLib::Point3d{
                std::array<double, 3>{{-wp[0], wp[1], -1.0}}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the pyramid.", face_id);
    }
}

MathLib::Point3d getBulkElementPoint(Tet const& /*tet*/,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp)
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
            return MathLib::Point3d{std::array<double, 3>{{0, wp[1], wp[0]}}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the tetrahedron.", face_id);
    }
}

}  // namespace MeshLib

static MathLib::Point3d getBulkElementPoint0D(MeshLib::Mesh const& mesh,
                                              std::size_t const bulk_element_id,
                                              std::size_t const bulk_face_id,
                                              MathLib::WeightedPoint const& wp)
{
    using namespace MeshLib;

    auto const* element = mesh.getElement(bulk_element_id);
    if (element->getCellType() == CellType::LINE2)
    {
        Line const& line(*dynamic_cast<Line const*>(element));
        return getBulkElementPoint(line, bulk_face_id, wp);
    }
    OGS_FATAL("Wrong cell type '{:s}' or functionality not yet implemented.",
              CellType2String(element->getCellType()));
}

static MathLib::Point3d getBulkElementPoint1D(MeshLib::Mesh const& mesh,
                                              std::size_t const bulk_element_id,
                                              std::size_t const bulk_face_id,
                                              MathLib::WeightedPoint const& wp)
{
    using namespace MeshLib;

    auto const* element = mesh.getElement(bulk_element_id);
    if (element->getCellType() == CellType::QUAD4)
    {
        Quad const& quad(*dynamic_cast<Quad const*>(element));
        return getBulkElementPoint(quad, bulk_face_id, wp);
    }
    if (element->getCellType() == CellType::TRI3)
    {
        Tri const& tri = *static_cast<Tri const*>(element);
        return getBulkElementPoint(tri, bulk_face_id, wp);
    }
    OGS_FATAL("Wrong cell type '{:s}' or functionality not yet implemented.",
              CellType2String(element->getCellType()));
}

static MathLib::Point3d getBulkElementPoint2D(MeshLib::Mesh const& mesh,
                                              std::size_t bulk_element_id,
                                              std::size_t bulk_face_id,
                                              MathLib::WeightedPoint const& wp)
{
    using namespace MeshLib;

    auto const* element = mesh.getElement(bulk_element_id);
    if (element->getCellType() == CellType::HEX8)
    {
        Hex const& hex = *static_cast<Hex const*>(element);
        return getBulkElementPoint(hex, bulk_face_id, wp);
    }
    if (element->getCellType() == CellType::PRISM6)
    {
        Prism const& prism = *static_cast<Prism const*>(element);
        return getBulkElementPoint(prism, bulk_face_id, wp);
    }
    if (element->getCellType() == CellType::PYRAMID5)
    {
        Pyramid const& pyramid = *static_cast<Pyramid const*>(element);
        return getBulkElementPoint(pyramid, bulk_face_id, wp);
    }
    if (element->getCellType() == CellType::TET4)
    {
        Tet const& tet = *static_cast<Tet const*>(element);
        return getBulkElementPoint(tet, bulk_face_id, wp);
    }
    OGS_FATAL("Wrong cell type '{:s}' or functionality not yet implemented.",
              CellType2String(element->getCellType()));
}

// TODO disable the 3d elements in the local assembler creator
static MathLib::Point3d getBulkElementPoint3D(
    MeshLib::Mesh const& /*mesh*/,
    std::size_t /*bulk_element_id*/,
    std::size_t /*bulk_face_id*/,
    MathLib::WeightedPoint const& /*wp*/)
{
    return MathLib::ORIGIN;
}

namespace MeshLib
{
MathLib::Point3d getBulkElementPoint(Mesh const& mesh,
                                     std::size_t const bulk_element_id,
                                     std::size_t const bulk_face_id,
                                     MathLib::WeightedPoint const& wp)
{
    switch (wp.getDimension())
    {
        case 0:
            return getBulkElementPoint0D(
                mesh, bulk_element_id, bulk_face_id, wp);
        case 1:
            return getBulkElementPoint1D(
                mesh, bulk_element_id, bulk_face_id, wp);
        case 2:
            return getBulkElementPoint2D(
                mesh, bulk_element_id, bulk_face_id, wp);
        case 3:
            return getBulkElementPoint3D(
                mesh, bulk_element_id, bulk_face_id, wp);
    }

    OGS_FATAL("Weighted point dimension of {} not supported.",
              wp.getDimension());
}
}  // namespace MeshLib
