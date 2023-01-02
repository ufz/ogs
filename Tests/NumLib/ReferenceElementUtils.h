/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/DynamicSpan.h"
#include "MeshLib/Elements/Element.h"

namespace ReferenceElementUtils
{

// Returns the coordinates as a span of dynamic size.
BaseLib::DynamicSpan<const std::array<double, 3>>
getNodeCoordsOfReferenceElement(MeshLib::CellType const cell_type);

std::shared_ptr<MeshLib::Element const> getReferenceElement(
    MeshLib::CellType const cell_type);

std::vector<std::array<double, 3>> getCoordsInReferenceElementForTest(
    MeshLib::Element const& element);

}  // namespace ReferenceElementUtils
