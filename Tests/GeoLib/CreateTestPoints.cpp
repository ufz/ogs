/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateTestPoints.h"

#include <map>
#include <memory>
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
                pnts.push_back(new GeoLib::Point(
                    i + shift[0], j + shift[1], k + shift[2], id));
                std::string const pnt_name(name + "-" + std::to_string(i) +
                                           "-" + std::to_string(j) + "-" +
                                           std::to_string(k));
                pnt_name_map.emplace(pnt_name, id);
            }
        }
    }

    geo_objs.addPointVec(std::move(pnts), name, std::move(pnt_name_map));
}
