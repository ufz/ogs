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

#include "MathLib/WeightedPoint.h"
#include "MeshLib/Elements/Elements.h"

namespace MeshLib
{

/// \page BulkMappingDocuPage
/// [Documentation](https://www.opengeosys.org/pdf/BulkElementMappings.pdf) for
/// the mapping of a point given in local coordinates of a boundary face/element
/// to the corresponding bulk element point.

/// Maps the given 0d end point \c wp of a line, to the corresponding 1d point
/// on the line. The result is given in local coordinates of the line element.
MathLib::Point3d getBulkElementPoint(Line const& line,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp);

/// Maps the given lower dimensional boundary point \c wp of a line, i.e. the 1d
/// integration point given in local coordinates of a line, to higher
/// dimensional point of the triangle face (defined by the triangle element and
/// the face id) also in local coordinates of the triangle element.
MathLib::Point3d getBulkElementPoint(Tri const& tri,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp);

/// Maps the given lower dimensional boundary point \c wp of a line, i.e. the 1d
/// integration point given in local coordinates of a line, to higher
/// dimensional point of the quad face (defined by the quad element and the face
/// id) also in local coordinates of the quad face.
MathLib::Point3d getBulkElementPoint(Quad const& quad,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp);

/// Maps the given 2d boundary point \c wp of a quad face, given in local
/// coordinates of the quad, to 3d point existing on a hexahedron face also in
/// local coordinates.
MathLib::Point3d getBulkElementPoint(Hex const& hex,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp);

/// Maps the given lower dimensional 2d boundary point \c wp of a triangle
/// element, given in local coordinates of the quad or triangle, to a 3d point
/// existing on a tet face also in local coordinates.
MathLib::Point3d getBulkElementPoint(Tet const& tet,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp);

/// Maps the given lower dimensional 2d boundary point \c wp of a quad or
/// triangle element, given in local coordinates of the quad or triangle, to a
/// 3d point existing on a prism face also in local coordinates.
MathLib::Point3d getBulkElementPoint(Prism const& prism,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp);

/// Maps the given lower dimensional 2d boundary point \c wp of a quad or
/// triangle element, given in local coordinates of the quad or triangle, to a
/// 3d point existing on a pyramid face also in local coordinates.
MathLib::Point3d getBulkElementPoint(Pyramid const& pyramid,
                                     std::size_t const face_id,
                                     MathLib::WeightedPoint const& wp);

/// Maps the given boundary point \c wp of a boundary element, given in local
/// coordinates of the boundary element, to 3d point existing on a bulk element
/// also in local coordinates.
MathLib::Point3d getBulkElementPoint(Mesh const& mesh,
                                     std::size_t const bulk_element_id,
                                     std::size_t const bulk_face_id,
                                     MathLib::WeightedPoint const& wp);
}  // namespace MeshLib
