/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateTestPoints.h"

#include <map>
#include <memory>
#include <random>
#include <vector>

void createSetOfTestPointsAndAssociatedNames(GeoLib::GEOObjects& geo_objs,
                                             std::string& name,
                                             std::size_t const pnts_per_edge,
                                             GeoLib::Point const& shift)
{
    std::vector<GeoLib::Point*> pnts;
    GeoLib::PointVec::NameIdMap pnt_name_map;

    for (std::size_t k(0); k < pnts_per_edge; k++)
    {
        const std::size_t k_offset(k * pnts_per_edge * pnts_per_edge);
        for (std::size_t j(0); j < pnts_per_edge; j++)
        {
            const std::size_t offset(j * pnts_per_edge + k_offset);
            for (std::size_t i(0); i < pnts_per_edge; i++)
            {
                std::size_t const id(i + offset);
                pnts.push_back(new GeoLib::Point(i + shift[0], j + shift[1],
                                                 k + shift[2], id));
                std::string const pnt_name(name + "-" + std::to_string(i) +
                                           "-" + std::to_string(j) + "-" +
                                           std::to_string(k));
                pnt_name_map.emplace(pnt_name, id);
            }
        }
    }

    geo_objs.addPointVec(std::move(pnts), name, std::move(pnt_name_map));
}

std::vector<GeoLib::Point*> createRandomPoints(
    std::size_t const number_of_random_points,
    std::array<double, 6> const& limits)
{
    static std::random_device rd;
    static std::mt19937 random_engine_mt19937(rd());
    std::normal_distribution<> normal_dist_x(limits[0], limits[1]);
    std::normal_distribution<> normal_dist_y(limits[2], limits[3]);
    std::normal_distribution<> normal_dist_z(limits[4], limits[5]);

    std::vector<GeoLib::Point*> random_points(number_of_random_points);

    auto random_point_generator = [&]()
    {
        return new GeoLib::Point(
            std::array{normal_dist_x(random_engine_mt19937),
                       normal_dist_y(random_engine_mt19937),
                       normal_dist_z(random_engine_mt19937)});
    };
    std::generate_n(random_points.begin(), number_of_random_points,
                    random_point_generator);

    return random_points;
}

std::vector<GeoLib::Point*> createPointsInGridArrangement(
    std::array<unsigned, 3> const& number_of_subdivisions,
    std::array<double, 3> const& distances,
    MathLib::Point3d const& origin)
{
    std::vector<GeoLib::Point*> points;
    points.reserve((number_of_subdivisions[0] + 1) *
                   (number_of_subdivisions[1] + 1) *
                   (number_of_subdivisions[2] + 1));

    for (std::size_t i = 0; i < number_of_subdivisions[2] + 1; i++)
    {
        const double z = origin[2] + distances[2] * i;
        for (std::size_t j = 0; j < number_of_subdivisions[1] + 1; j++)
        {
            const double y = origin[1] + distances[1] * j;
            for (std::size_t k = 0; k < number_of_subdivisions[0] + 1; k++)
            {
                points.push_back(
                    new GeoLib::Point(origin[0] + distances[0] * k, y, z));
            }
        }
    }
    return points;
}
