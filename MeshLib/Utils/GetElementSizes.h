/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 2, 2022, 3:05 PM
 */

#pragma once

#include <vector>

namespace MeshLib
{
class Mesh;
/// Returns the maximum lengths of the edges for each element.
std::vector<double> getElementSizes(Mesh const& mesh);
}  // namespace MeshLib
