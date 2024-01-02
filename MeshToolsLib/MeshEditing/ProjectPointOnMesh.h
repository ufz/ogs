/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

namespace MeshToolsLib
{
namespace ProjectPointOnMesh
{

/// Returns the element in which the given node is located when projected onto a
/// mesh, or nullptr if no such element was found.
MeshLib::Element const* getProjectedElement(
    std::vector<const MeshLib::Element*> const& elements,
    MathLib::Point3d const& node);

/// Returns the z-coordinate of a point projected onto the plane defined by a
/// mesh element.
double getElevation(MeshLib::Element const& element,
                    MathLib::Point3d const& node);

}  // namespace ProjectPointOnMesh

}  // namespace MeshToolsLib
