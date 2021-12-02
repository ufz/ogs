/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "GeoLib/Point.h"

namespace GeoLib
{
/// Function generates equidistant points in the interval [begin, end]
/// according to the given number of subdivisions.
/// In case of zero subdivisions a vector containing the begin and the end point
/// is returned.
std::vector<GeoLib::Point*> generateEquidistantPoints(
    MathLib::Point3d const& begin, MathLib::Point3d const& end,
    int const number_of_subdivisions);
}  // namespace GeoLib
