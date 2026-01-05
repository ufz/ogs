// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <span>

#include "MeshLib/Elements/Element.h"
#include "NumLib/Fem/CoordinatesMapping/NaturalNodeCoordinates.h"

namespace ReferenceElementUtils
{

// Returns the coordinates as a span of known extent.
template <typename MeshElementType>
auto getNodeCoordsOfReferenceElement()
{
    return std::span{NumLib::NaturalCoordinates<MeshElementType>::coordinates};
}

// Returns the coordinates as a span of dynamic size.
std::span<const std::array<double, 3>> getNodeCoordsOfReferenceElement(
    MeshLib::CellType const cell_type);

std::shared_ptr<MeshLib::Element const> getReferenceElement(
    MeshLib::CellType const cell_type);

std::vector<std::array<double, 3>> getCoordsInReferenceElementForTest(
    MeshLib::Element const& element);

}  // namespace ReferenceElementUtils
