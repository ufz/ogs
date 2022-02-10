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

/// Maps the given boundary point \c wp of a boundary element, given in local
/// coordinates of the boundary element, to 3d point existing on a bulk element
/// also in local coordinates.
MathLib::Point3d getBulkElementPoint(Mesh const& mesh,
                                     std::size_t const bulk_element_id,
                                     std::size_t const bulk_face_id,
                                     MathLib::WeightedPoint const& wp);
}  // namespace MeshLib
