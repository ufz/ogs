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
    : mesh_(mesh), mshNodeSearcher_(mshNodeSearcher)
{}

BoundaryElementsSearcher::~BoundaryElementsSearcher()
{
    for (auto p : boundary_elements_at_point_)
    {
        delete p;
    }
    for (auto p : boundary_elements_along_polylines_)
    {
        delete p;
    }
    for (auto p : boundary_elements_along_surfaces_)
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
    for (auto const& boundaryElements : boundary_elements_at_point_)
    {
        if (boundaryElements->getPoint() == point)
        {
            return boundaryElements->getBoundaryElements();
        }
    }

    // create new boundary elements at points.
    boundary_elements_at_point_.push_back(new BoundaryElementsAtPoint(
        mesh_, mshNodeSearcher_, point, multiple_nodes_allowed));
    return boundary_elements_at_point_.back()->getBoundaryElements();
}

std::vector<MeshLib::Element*> const&
BoundaryElementsSearcher::getBoundaryElementsAlongPolyline(
    GeoLib::Polyline const& polyline)
{
    // look for already saved polylines and return if found.
    for (auto const& boundary_elements : boundary_elements_along_polylines_)
    {
        if (&boundary_elements->getPolyline() == &polyline)
        {
            return boundary_elements->getBoundaryElements();
        }
    }

    // create new boundary elements at points.
    boundary_elements_along_polylines_.push_back(
        new BoundaryElementsAlongPolyline(mesh_, mshNodeSearcher_, polyline));
    return boundary_elements_along_polylines_.back()->getBoundaryElements();
}

std::vector<MeshLib::Element*> const&
BoundaryElementsSearcher::getBoundaryElementsOnSurface(
    GeoLib::Surface const& surface)
{
    // look for already saved surfaces and return if found.
    for (auto const& boundary_elements : boundary_elements_along_surfaces_)
    {
        if (&boundary_elements->getSurface() == &surface)
        {
            return boundary_elements->getBoundaryElements();
        }
    }

    boundary_elements_along_surfaces_.push_back(
            new BoundaryElementsOnSurface(mesh_, mshNodeSearcher_, surface));
    return boundary_elements_along_surfaces_.back()->getBoundaryElements();
}

} // end namespace MeshGeoToolsLib

