// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
