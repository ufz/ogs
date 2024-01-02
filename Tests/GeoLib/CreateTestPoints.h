/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <array>
#include <string>

#include "GeoLib/GEOObjects.h"

void createSetOfTestPointsAndAssociatedNames(GeoLib::GEOObjects& geo_objs,
                                             std::string& name,
                                             std::size_t const pnts_per_edge,
                                             GeoLib::Point const& shift);

std::vector<GeoLib::Point*> createRandomPoints(
    std::size_t const number_of_random_points,
    std::array<double, 6> const& limits);

std::vector<GeoLib::Point*> createPointsInGridArrangement(
    std::array<unsigned, 3> const& number_of_subdivisions,
    std::array<double, 3> const& distances,
    MathLib::Point3d const& origin);
