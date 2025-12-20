// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
