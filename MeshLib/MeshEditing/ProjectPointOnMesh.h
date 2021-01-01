/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
 MeshLib::Element const* getProjectedElement(
    std::vector<const MeshLib::Element*> const& elements,
    MeshLib::Node const& node);

/// Returns the z-coordinate of a point projected onto the plane defined by a
/// mesh element.
double getElevation(MeshLib::Element const& element,
                    MeshLib::Node const& node);

}  // namespace ProjectPointOnMesh

}  // end namespace MeshLib
