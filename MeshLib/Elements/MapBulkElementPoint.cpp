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

namespace
{
/// Maps the given 0d end point \c wp of a line to the corresponding 1d point
/// on the line.
///
/// The result is given in local coordinates of the line element.
MathLib::Point3d getBulkElementPointLine(std::size_t const face_id)
{
    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{{-1, 0, 0}};
        case 1:
            return MathLib::Point3d{{1, 0, 0}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for a line.", face_id);
    }
}

/// Maps the given lower dimensional boundary point \c wp of a line, i.e. the 1d
/// integration point given in local coordinates of a line, to higher
/// dimensional point of the triangle face (defined by the triangle element and
/// the face id) also in local coordinates of the triangle element.
///
/// \param face_id the id of the triangle face the point will be mapped on
/// \param wp the integration point of the lower dimensional element
/// \return the mapped point
MathLib::Point3d getBulkElementPointTri(std::size_t const face_id,
                                        MathLib::WeightedPoint const& wp)
{
    // wp is a local coordinate of the 1D line face of the triangle. I.e., wp
    // has only one relevant coordinate. This coordinate runs from -1 to 1.
    auto const r = 0.5 * (wp[0] + 1);  // in [0, 1]

    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{{r, 0, 0}};
        case 1:
            return MathLib::Point3d{{1 - r, r, 0}};
        case 2:
            return MathLib::Point3d{{0, 1 - r, 0}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the tri.", face_id);
    }
}

/// Maps the given lower dimensional boundary point \c wp of a line, i.e. the 1d
/// integration point given in local coordinates of a line, to higher
/// dimensional point of the quad face (defined by the quad element and the face
/// id) also in local coordinates of the quad face.
///
/// \param face_id the id of the quad face the point will be mapped on
/// \param wp the integration point of the lower dimensional element
/// \return the mapped point
MathLib::Point3d getBulkElementPointQuad(std::size_t const face_id,
                                         MathLib::WeightedPoint const& wp)
{
    auto const r = wp[0];

    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{{-r, +1, 0}};
        case 1:
            return MathLib::Point3d{{-1, -r, 0}};
        case 2:
            return MathLib::Point3d{{+r, -1, 0}};
        case 3:
            return MathLib::Point3d{{+1, +r, 0}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the quad.", face_id);
    }
}

/// Maps the given 2d boundary point \c wp of a hexahedron face to the
/// corresponding 3d point of the hexahedron.
///
/// The input and output coordinates are natural coordinates of the quad and
/// hex, respectively.
MathLib::Point3d getBulkElementPointHex(std::size_t const face_id,
                                        MathLib::WeightedPoint const& wp)
{
    auto const r = wp[0];
    auto const s = wp[1];

    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{{-s, -r, -1}};
        case 1:
            return MathLib::Point3d{{-r, -1, -s}};
        case 2:
            return MathLib::Point3d{{+1, -r, -s}};
        case 3:
            return MathLib::Point3d{{+r, +1, -s}};
        case 4:
            return MathLib::Point3d{{-1, +r, -s}};
        case 5:
            return MathLib::Point3d{{-r, -s, +1}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the hexahedron.", face_id);
    }
}

/// Maps the given lower dimensional 2d boundary point \c wp of a quad or
/// triangle element, given in local coordinates of the quad or triangle, to a
/// 3d point existing on a prism face also in local coordinates.
MathLib::Point3d getBulkElementPointPrism(std::size_t const face_id,
                                          MathLib::WeightedPoint const& wp)
{
    auto const r = wp[0];
    auto const s = wp[1];

    switch (face_id)
    {
        case 0:
            // tri, top; r, s in [0, 1]
            return MathLib::Point3d{{s, r, -1}};
        case 1:
        {
            // quad; r, s in [-1, 1]
            auto const r01 = 0.5 * (r + 1);  // in [0, 1]
            return MathLib::Point3d{{1 - r01, 0, -s}};
        }
        case 2:
        {
            // quad; r, s in [-1, 1]
            auto const r01 = 0.5 * (r + 1);  // in [0, 1]
            return MathLib::Point3d{{r01, 1 - r01, -s}};
        }
        case 3:
        {
            // quad; r, s in [-1, 1]
            auto const r01 = 0.5 * (r + 1);  // in [0, 1]
            return MathLib::Point3d{{0, r01, -s}};
        }
        case 4:
            // tri, bottom; r, s in [0, 1]
            return MathLib::Point3d{{r, s, 1}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the prism.", face_id);
    }
}

/// Maps the given lower dimensional 2d boundary point \c wp of a quad or
/// triangle element, given in local coordinates of the quad or triangle, to a
/// 3d point existing on a pyramid face also in local coordinates.
MathLib::Point3d getBulkElementPointPyramid(std::size_t const face_id,
                                            MathLib::WeightedPoint const& wp)
{
    auto const r = wp[0];  // in [0, 1] on tri faces, in [-1, 1] on quad face
    auto const s = wp[1];  // in [0, 1] on tri faces, in [-1, 1] on quad face

    // for tri phases, note: there r + s == 1, i.e., s == 1 => r == 0
    auto const z = 2 * s - 1;  // in [-1, 1] on tri faces
    auto const w = s == 1 ? -1 : 2 * r / (1 - s) - 1;

    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{{+w, -1.0, z}};
        case 1:
            return MathLib::Point3d{{+1.0, +w, z}};
        case 2:
            return MathLib::Point3d{{-w, +1.0, z}};
        case 3:
            return MathLib::Point3d{{-1.0, -w, z}};
        case 4:
            return MathLib::Point3d{{-s, -r, -1.0}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the pyramid.", face_id);
    }
}

/// Maps the given lower dimensional 2d boundary point \c wp of a triangle
/// element, given in local coordinates of the quad or triangle, to a 3d point
/// existing on a tet face also in local coordinates.
MathLib::Point3d getBulkElementPointTet(std::size_t const face_id,
                                        MathLib::WeightedPoint const& wp)
{
    auto const r = wp[0];
    auto const s = wp[1];
    auto const t = 1 - r - s;

    switch (face_id)
    {
        case 0:
            return MathLib::Point3d{{s, r, 0}};
        case 1:
            return MathLib::Point3d{{r, 0, s}};
        case 2:
            return MathLib::Point3d{{t, r, s}};
        case 3:
            return MathLib::Point3d{{0, t, s}};
        default:
            OGS_FATAL("Invalid face id '{:d}' for the tetrahedron.", face_id);
    }
}

}  // namespace

namespace MeshLib
{
MathLib::Point3d getBulkElementPoint(
    MeshLib::CellType const bulk_element_cell_type,
    std::size_t const bulk_face_id,
    MathLib::WeightedPoint const& point_on_face)
{
    if (point_on_face.getDimension() == 3)
    {
        // TODO disable the 3d elements in the local assembler creator
        return MathLib::ORIGIN;
    }

    switch (bulk_element_cell_type)
    {
            // 3D bulk elements
        case CellType::HEX8:
            return getBulkElementPointHex(bulk_face_id, point_on_face);
        case CellType::PRISM6:
            return getBulkElementPointPrism(bulk_face_id, point_on_face);
        case CellType::PYRAMID5:
            return getBulkElementPointPyramid(bulk_face_id, point_on_face);
        case CellType::TET4:
            return getBulkElementPointTet(bulk_face_id, point_on_face);
            // 2D bulk elements
        case CellType::QUAD4:
            return getBulkElementPointQuad(bulk_face_id, point_on_face);
        case CellType::TRI3:
            return getBulkElementPointTri(bulk_face_id, point_on_face);
            // 1D bulk elements
        case CellType::LINE2:
            return getBulkElementPointLine(bulk_face_id);
        default:
            OGS_FATAL(
                "Wrong cell type '{:s}' or functionality not yet implemented.",
                CellType2String(bulk_element_cell_type));
    }
}
}  // namespace MeshLib
