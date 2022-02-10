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

/// Maps the given \c point_on_face on the face of the given \c bulk_element to
/// the given point in the \c bulk_element.
///
/// The input and output coordinates are natural coordinates of the surface and
/// bulk element, respectively. I.e., the output point has one coordinate more
/// than the input point.
MathLib::Point3d getBulkElementPoint(
    MeshLib::Element const& bulk_element,
    std::size_t const bulk_face_id,
    MathLib::WeightedPoint const& point_on_face);
}  // namespace MeshLib
