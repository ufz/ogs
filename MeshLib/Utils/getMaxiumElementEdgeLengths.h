// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

namespace MeshLib
{
class Element;
/// Returns the maximum lengths of the edges for each element.
std::vector<double> getMaxiumElementEdgeLengths(
    std::vector<Element*> const& elements);
}  // namespace MeshLib
