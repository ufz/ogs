/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "BoundaryElementsSearcher.h"

#include "GeoLib/GeoObject.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Point.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/BoundaryElementsAtPoint.h"
#include "MeshGeoToolsLib/BoundaryElementsAlongPolyline.h"
#include "MeshGeoToolsLib/BoundaryElementsOnSurface.h"


namespace MeshGeoToolsLib
{
BoundaryElementsSearcher::BoundaryElementsSearcher(
    MeshLib::Mesh const& mesh, MeshNodeSearcher const& mshNodeSearcher)
    : _mesh(mesh), _mshNodeSearcher(mshNodeSearcher)
{}

BoundaryElementsSearcher::~BoundaryElementsSearcher()
{
    for (auto p : _boundary_elements_at_point)
        delete p;
    for (auto p : _boundary_elements_along_polylines)
        delete p;
    for (auto p : _boundary_elements_along_surfaces)
        delete p;
}

std::vector<MeshLib::Element*> const& BoundaryElementsSearcher::getBoundaryElements(GeoLib::GeoObject const& geoObj)
{
    switch (geoObj.getGeoType()) {
    case GeoLib::GEOTYPE::POINT:
        return this->getBoundaryElementsAtPoint(*dynamic_cast<const GeoLib::Point*>(&geoObj));
        break;
    case GeoLib::GEOTYPE::POLYLINE:
        return this->getBoundaryElementsAlongPolyline(*dynamic_cast<const GeoLib::Polyline*>(&geoObj));
        break;
    case GeoLib::GEOTYPE::SURFACE:
        return this->getBoundaryElementsOnSurface(*dynamic_cast<const GeoLib::Surface*>(&geoObj));
        break;
    default:
        const static std::vector<MeshLib::Element*> dummy(0);
        return dummy;
    }
}

std::vector<MeshLib::Element*> const&
BoundaryElementsSearcher::getBoundaryElementsAtPoint(GeoLib::Point const& point)
{
    // look for already saved points and return if found.
    for (auto const& boundaryElements : _boundary_elements_at_point)
    {
        if (boundaryElements->getPoint() == point)
            return boundaryElements->getBoundaryElements();
    }

    // create new boundary elements at points.
    _boundary_elements_at_point.push_back(
        new BoundaryElementsAtPoint(_mesh, _mshNodeSearcher, point));
    return _boundary_elements_at_point.back()->getBoundaryElements();
}

std::vector<MeshLib::Element*> const& BoundaryElementsSearcher::getBoundaryElementsAlongPolyline(GeoLib::Polyline const& ply)
{
    std::vector<BoundaryElementsAlongPolyline*>::const_iterator it(_boundary_elements_along_polylines.begin());
    for (; it != _boundary_elements_along_polylines.end(); ++it) {
        if (&(*it)->getPolyline() == &ply) {
            // we calculated mesh nodes for this polyline already
            return (*it)->getBoundaryElements();
        }
    }

    _boundary_elements_along_polylines.push_back(
            new BoundaryElementsAlongPolyline(_mesh, _mshNodeSearcher, ply));
    return _boundary_elements_along_polylines.back()->getBoundaryElements();
}

std::vector<MeshLib::Element*> const& BoundaryElementsSearcher::getBoundaryElementsOnSurface(GeoLib::Surface const& sfc)
{
    std::vector<BoundaryElementsOnSurface*>::const_iterator it(_boundary_elements_along_surfaces.begin());
    for (; it != _boundary_elements_along_surfaces.end(); ++it) {
        if (&(*it)->getSurface() == &sfc) {
            // we calculated mesh nodes for this surface already
            return (*it)->getBoundaryElements();
        }
    }

    _boundary_elements_along_surfaces.push_back(
            new BoundaryElementsOnSurface(_mesh, _mshNodeSearcher, sfc));
    return _boundary_elements_along_surfaces.back()->getBoundaryElements();
}

} // end namespace MeshGeoToolsLib

