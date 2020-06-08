/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
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
    if (auto const it = find_if(cbegin(_boundary_elements_at_point),
                                cend(_boundary_elements_at_point),
                                [&](auto const& boundary_elements) {
                                    return boundary_elements->getPoint() ==
                                           point;
                                });
        it != cend(_boundary_elements_at_point))
    {
        return (*it)->getBoundaryElements();
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
    if (auto const it = find_if(cbegin(_boundary_elements_along_polylines),
                                cend(_boundary_elements_along_polylines),
                                [&](auto const& boundary_elements) {
                                    return &boundary_elements->getPolyline() ==
                                           &polyline;
                                });
        it != cend(_boundary_elements_along_polylines))
    {
        return (*it)->getBoundaryElements();
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
    if (auto const it = find_if(cbegin(_boundary_elements_along_surfaces),
                                cend(_boundary_elements_along_surfaces),
                                [&](auto const& boundary_elements) {
                                    return &boundary_elements->getSurface() ==
                                           &surface;
                                });
        it != cend(_boundary_elements_along_surfaces))
    {
        return (*it)->getBoundaryElements();
    }

    _boundary_elements_along_surfaces.push_back(
            new BoundaryElementsOnSurface(_mesh, _mshNodeSearcher, surface));
    return _boundary_elements_along_surfaces.back()->getBoundaryElements();
}

} // end namespace MeshGeoToolsLib

