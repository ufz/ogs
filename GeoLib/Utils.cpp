/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Utils.h"

#include "BaseLib/Error.h"

namespace GeoLib
{
std::vector<GeoLib::Point*> generateEquidistantPoints(
    MathLib::Point3d const& begin, MathLib::Point3d const& end,
    int const number_of_subdivisions)
{
    if (number_of_subdivisions < 0)
    {
        OGS_FATAL(
            "generateEquidistantPoints: number of subdivisions is required to "
            "be non-negative.");
    }

    auto const start = Eigen::Map<Eigen::Vector3d const>(begin.data());
    auto const stop = Eigen::Map<Eigen::Vector3d const>(end.data());
    auto const delta = (stop - start) / (number_of_subdivisions + 1);

    std::vector<GeoLib::Point*> points;

    for (int i = 0; i <= number_of_subdivisions; ++i)
    {
        auto const p = start + i * delta;
        points.push_back(new GeoLib::Point{p[0], p[1], p[2]});
    }
    points.push_back(new GeoLib::Point{stop[0], stop[1], stop[2]});

    return points;
}

}  // namespace GeoLib
