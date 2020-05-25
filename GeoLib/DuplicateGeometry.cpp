/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DuplicateGeometry.h"

#include <map>
#include <utility>
#include "BaseLib/Logging.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"
#include "GeoLib/Triangle.h"


namespace GeoLib
{
DuplicateGeometry::DuplicateGeometry(GeoLib::GEOObjects& geo_objects,
                                     std::string const& input_name,
                                     std::string output_name)
    : output_name_(std::move(output_name)), geo_objects_(geo_objects)
{
    duplicate(input_name);
}

void DuplicateGeometry::duplicate(std::string const& input_name)
{
    std::vector<GeoLib::Point*> const*const pnts (geo_objects_.getPointVec(input_name));
    if (pnts == nullptr)
    {
        ERR("Geometry '{:s}' not found.", input_name);
        return;
    }

    auto new_pnts = std::make_unique<std::vector<GeoLib::Point*>>();
    new_pnts->reserve(pnts->size());
    std::transform(pnts->cbegin(), pnts->cend(), std::back_inserter(*new_pnts),
        [](GeoLib::Point* point) { return new GeoLib::Point(*point); });
    auto pnt_name_id_map = std::make_unique<std::map<std::string, std::size_t>>(
        geo_objects_.getPointVecObj(input_name)->getNameIDMapBegin(),
        geo_objects_.getPointVecObj(input_name)->getNameIDMapEnd());
    geo_objects_.addPointVec(std::move(new_pnts), output_name_,
                             std::move(pnt_name_id_map));

    std::vector<GeoLib::Polyline*> const* plys (geo_objects_.getPolylineVec(input_name));
    if (plys)
    {
        auto new_plys = copyPolylinesVector(*plys);
        auto ply_name_id_map =
            std::make_unique<std::map<std::string, std::size_t>>(
                geo_objects_.getPolylineVecObj(input_name)->getNameIDMapBegin(),
                geo_objects_.getPolylineVecObj(input_name)->getNameIDMapEnd());
        geo_objects_.addPolylineVec(std::move(new_plys), output_name_,
                                    std::move(ply_name_id_map));
    }

    std::vector<GeoLib::Surface*> const* sfcs (geo_objects_.getSurfaceVec(input_name));
    if (sfcs)
    {
        auto new_sfcs = copySurfacesVector(*sfcs);
        auto sfc_name_id_map =
            std::make_unique<std::map<std::string, std::size_t>>(
                geo_objects_.getSurfaceVecObj(input_name)->getNameIDMapBegin(),
                geo_objects_.getSurfaceVecObj(input_name)->getNameIDMapEnd());
        geo_objects_.addSurfaceVec(std::move(new_sfcs), output_name_,
                                   std::move(sfc_name_id_map));
    }
}

std::unique_ptr<std::vector<GeoLib::Polyline*>> DuplicateGeometry::copyPolylinesVector(
    std::vector<GeoLib::Polyline*> const& polylines) const
{
    std::size_t const n_plys = polylines.size();
    auto new_lines =
        std::make_unique<std::vector<GeoLib::Polyline*>>(n_plys, nullptr);

    for (std::size_t i=0; i<n_plys; ++i)
    {
        if (polylines[i] == nullptr)
        {
            continue;
        }
        (*new_lines)[i] = new GeoLib::Polyline(*geo_objects_.getPointVec(output_name_));
        std::size_t const nLinePnts (polylines[i]->getNumberOfPoints());
        for (std::size_t j = 0; j < nLinePnts; ++j)
        {
            (*new_lines)[i]->addPoint(polylines[i]->getPointID(j));
        }
    }
    return new_lines;
}

std::unique_ptr<std::vector<Surface*>> DuplicateGeometry::copySurfacesVector(
    std::vector<Surface*> const& surfaces) const
{
    std::size_t const n_sfc = surfaces.size();
    auto new_surfaces =
        std::make_unique<std::vector<GeoLib::Surface*>>(n_sfc, nullptr);

    for (std::size_t i=0; i<n_sfc; ++i)
    {
        if (surfaces[i] == nullptr)
        {
            continue;
        }
        (*new_surfaces)[i] = new GeoLib::Surface(*geo_objects_.getPointVec(output_name_));

        std::size_t const n_tris (surfaces[i]->getNumberOfTriangles());
        for (std::size_t j=0; j<n_tris; ++j)
        {
            GeoLib::Triangle const* t = (*surfaces[i])[j];
            (*new_surfaces)[i]->addTriangle(
                t->getPoint(0)->getID(), t->getPoint(1)->getID(), t->getPoint(2)->getID());
        }
    }
    return new_surfaces;
}

std::vector<GeoLib::Point*>& DuplicateGeometry::getPointVectorCopy()
{
    return const_cast<std::vector<GeoLib::Point*>&>(*geo_objects_.getPointVec(output_name_));
}

std::vector<GeoLib::Polyline*>& DuplicateGeometry::getPolylineVectorCopy()
{
    return const_cast<std::vector<GeoLib::Polyline*>&>(*geo_objects_.getPolylineVec(output_name_));
}

std::vector<GeoLib::Surface*>& DuplicateGeometry::getSurfaceVectorCopy()
{
    return const_cast<std::vector<GeoLib::Surface*>&>(*geo_objects_.getSurfaceVec(output_name_));
}

}  // namespace GeoLib
