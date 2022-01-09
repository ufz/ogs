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

#include <vector>
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib
{
namespace ProjectPointOnMesh
{

/// Returns the element in which the given node is located when projected onto a
/// mesh, or nullptr if no such element was found.
Element const* getProjectedElement(std::vector<const Element*> const& elements,
                                   Node const& node);

/// Returns the z-coordinate of a point projected onto the plane defined by a
/// mesh element.
double getElevation(Element const& element, Node const& node);

}  // namespace ProjectPointOnMesh

}  // end namespace MeshLib
