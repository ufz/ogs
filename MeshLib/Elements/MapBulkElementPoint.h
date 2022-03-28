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
#include "MeshLib/Elements/Element.h"

namespace MeshLib
{
/// \page BulkMappingDocuPage
/// [Documentation](https://www.opengeosys.org/pdf/BulkElementMappings.pdf) for
/// the mapping of a point given in local coordinates of a boundary face/element
/// to the corresponding bulk element point.

/// Maps the given \c point_on_face on the face a bulk element of type \c
/// bulk_element_cell_type to the given point in the bulk element.
///
/// The input and output coordinates are natural coordinates of the surface and
/// bulk element, respectively. I.e., the output point has one coordinate more
/// than the input point.
///
MathLib::Point3d getBulkElementPoint(
    MeshLib::CellType const bulk_element_cell_type,
    std::size_t const bulk_face_id,
    MathLib::WeightedPoint const& point_on_face);

/// Overload provided for convenience.
inline MathLib::Point3d getBulkElementPoint(
    MeshLib::Element const& bulk_element,
    std::size_t const bulk_face_id,
    MathLib::WeightedPoint const& point_on_face)
{
    return getBulkElementPoint(
        bulk_element.getCellType(), bulk_face_id, point_on_face);
}
}  // namespace MeshLib
