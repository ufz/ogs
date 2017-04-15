/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DuplicateGeometry.h"

#include <logog/include/logog.hpp>
#include <utility>

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
    : _output_name(std::move(output_name)), _geo_objects(geo_objects)
{
    duplicate(input_name);
}

void DuplicateGeometry::duplicate(std::string const& input_name)
{
    std::vector<GeoLib::Point*> const*const pnts (_geo_objects.getPointVec(input_name));
    if (pnts == nullptr)
    {
        ERR("Geometry \"%s\" not found.", input_name.c_str());
        return;
    }

    std::unique_ptr< std::vector<GeoLib::Point*> > new_pnts (new std::vector<GeoLib::Point*>);
    new_pnts->reserve(pnts->size());
    std::transform(pnts->cbegin(), pnts->cend(), std::back_inserter(*new_pnts),
        [](GeoLib::Point* point) { return new GeoLib::Point(*point); });
    _geo_objects.addPointVec(std::move(new_pnts), _output_name);

    std::vector<GeoLib::Polyline*> const* plys (_geo_objects.getPolylineVec(input_name));
    if (plys)
    {
        auto new_plys = copyPolylinesVector(*plys);
        _geo_objects.addPolylineVec(std::move(new_plys), _output_name);
    }

    std::vector<GeoLib::Surface*> const* sfcs (_geo_objects.getSurfaceVec(input_name));
    if (sfcs)
    {
        auto new_sfcs = copySurfacesVector(*sfcs);
        _geo_objects.addSurfaceVec(std::move(new_sfcs), _output_name);
    }
}

std::unique_ptr<std::vector<GeoLib::Polyline*>> DuplicateGeometry::copyPolylinesVector(
    std::vector<GeoLib::Polyline*> const& polylines) const
{
    std::size_t const n_plys = polylines.size();
    std::unique_ptr<std::vector<GeoLib::Polyline*>> new_lines
        (new std::vector<Polyline*>(n_plys, nullptr));

    for (std::size_t i=0; i<n_plys; ++i)
    {
        if (polylines[i] == nullptr)
            continue;
        (*new_lines)[i] = new GeoLib::Polyline(*_geo_objects.getPointVec(_output_name));
        std::size_t const nLinePnts (polylines[i]->getNumberOfPoints());
        for (std::size_t j=0; j<nLinePnts; ++j)
            (*new_lines)[i]->addPoint(polylines[i]->getPointID(j));
    }
    return new_lines;
}

std::unique_ptr<std::vector<Surface*>> DuplicateGeometry::copySurfacesVector(
    std::vector<Surface*> const& surfaces) const
{
    std::size_t const n_sfc = surfaces.size();
    std::unique_ptr<std::vector<GeoLib::Surface*>> new_surfaces
        (new std::vector<Surface*>(n_sfc, nullptr));

    for (std::size_t i=0; i<n_sfc; ++i)
    {
        if (surfaces[i] == nullptr)
            continue;
        (*new_surfaces)[i] = new GeoLib::Surface(*_geo_objects.getPointVec(_output_name));

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
    return const_cast<std::vector<GeoLib::Point*>&>(*_geo_objects.getPointVec(_output_name));
}

std::vector<GeoLib::Polyline*>& DuplicateGeometry::getPolylineVectorCopy()
{
    return const_cast<std::vector<GeoLib::Polyline*>&>(*_geo_objects.getPolylineVec(_output_name));
}

std::vector<GeoLib::Surface*>& DuplicateGeometry::getSurfaceVectorCopy()
{
    return const_cast<std::vector<GeoLib::Surface*>&>(*_geo_objects.getSurfaceVec(_output_name));
}

}
