/**
 * @copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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
    {
        delete p;
    }
    for (auto p : _boundary_elements_along_polylines)
    {
        delete p;
    }
    for (auto p : _boundary_elements_along_surfaces)
    {
        delete p;
    }
}

std::vector<MeshLib::Element*> const&
BoundaryElementsSearcher::getBoundaryElements(GeoLib::GeoObject const& geoObj,
                                              bool const multiple_nodes_allowed)
{
    switch (geoObj.getGeoType()) {
    case GeoLib::GEOTYPE::POINT:
        return this->getBoundaryElementsAtPoint(
            *dynamic_cast<const GeoLib::Point*>(&geoObj),
            multiple_nodes_allowed);
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
BoundaryElementsSearcher::getBoundaryElementsAtPoint(
    GeoLib::Point const& point, bool const multiple_nodes_allowed)
{
    // look for already saved points and return if found.
    for (auto const& boundaryElements : _boundary_elements_at_point)
    {
        if (boundaryElements->getPoint() == point)
        {
            return boundaryElements->getBoundaryElements();
        }
    }

    // create new boundary elements at points.
    _boundary_elements_at_point.push_back(new BoundaryElementsAtPoint(
        _mesh, _mshNodeSearcher, point, multiple_nodes_allowed));
    return _boundary_elements_at_point.back()->getBoundaryElements();
}

std::vector<MeshLib::Element*> const&
BoundaryElementsSearcher::getBoundaryElementsAlongPolyline(
    GeoLib::Polyline const& polyline)
{
    // look for already saved polylines and return if found.
    for (auto const& boundary_elements : _boundary_elements_along_polylines)
    {
        if (&boundary_elements->getPolyline() == &polyline)
        {
            return boundary_elements->getBoundaryElements();
        }
    }

    // create new boundary elements at points.
    _boundary_elements_along_polylines.push_back(
        new BoundaryElementsAlongPolyline(_mesh, _mshNodeSearcher, polyline));
    return _boundary_elements_along_polylines.back()->getBoundaryElements();
}

std::vector<MeshLib::Element*> const&
BoundaryElementsSearcher::getBoundaryElementsOnSurface(
    GeoLib::Surface const& surface)
{
    // look for already saved surfaces and return if found.
    for (auto const& boundary_elements : _boundary_elements_along_surfaces)
    {
        if (&boundary_elements->getSurface() == &surface)
        {
            return boundary_elements->getBoundaryElements();
        }
    }

    _boundary_elements_along_surfaces.push_back(
            new BoundaryElementsOnSurface(_mesh, _mshNodeSearcher, surface));
    return _boundary_elements_along_surfaces.back()->getBoundaryElements();
}

} // end namespace MeshGeoToolsLib

