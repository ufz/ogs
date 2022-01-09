/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BoundaryElementsSearcher.h"

#include "GeoLib/GeoObject.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"
#include "MeshGeoToolsLib/BoundaryElementsAlongPolyline.h"
#include "MeshGeoToolsLib/BoundaryElementsAtPoint.h"
#include "MeshGeoToolsLib/BoundaryElementsOnSurface.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Point.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshGeoToolsLib
{
BoundaryElementsSearcher::BoundaryElementsSearcher(
    MeshLib::Mesh const& mesh, MeshNodeSearcher const& mshNodeSearcher)
    : _mesh(mesh), _mshNodeSearcher(mshNodeSearcher)
{
}

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

template <typename CacheType, typename GeometryType>
std::vector<MeshLib::Element*> const& getBoundaryElements(
    std::vector<CacheType*>& cached_elements,
    std::function<GeometryType(CacheType const&)> getCachedItem,
    GeometryType const& item, MeshLib::Mesh const& mesh,
    MeshNodeSearcher const& mesh_node_searcher,
    [[maybe_unused]] bool const multiple_nodes_allowed)
{
    if (auto const it = find_if(cbegin(cached_elements), cend(cached_elements),
                                [&](auto const& element)
                                { return getCachedItem(*element) == item; });
        it != cend(cached_elements))
    {
        return (*it)->getBoundaryElements();
    }
    // create new boundary elements
    if constexpr (std::is_convertible<GeometryType, GeoLib::Point>::value)
    {
        cached_elements.push_back(new CacheType(mesh, mesh_node_searcher, item,
                                                multiple_nodes_allowed));
    }
    else
    {
        cached_elements.push_back(
            new CacheType(mesh, mesh_node_searcher, item));
    }
    return cached_elements.back()->getBoundaryElements();
}

std::vector<MeshLib::Element*> const&
BoundaryElementsSearcher::getBoundaryElements(GeoLib::GeoObject const& geoObj,
                                              bool const multiple_nodes_allowed)
{
    switch (geoObj.getGeoType())
    {
        case GeoLib::GEOTYPE::POINT:
        {
            std::function<GeoLib::Point(BoundaryElementsAtPoint const&)>
                get_cached_item = &BoundaryElementsAtPoint::getPoint;
            return MeshGeoToolsLib::getBoundaryElements(
                _boundary_elements_at_point, get_cached_item,
                *dynamic_cast<const GeoLib::Point*>(&geoObj), _mesh,
                _mshNodeSearcher, multiple_nodes_allowed);
        }
        break;
        case GeoLib::GEOTYPE::POLYLINE:
        {
            std::function<GeoLib::Polyline(
                BoundaryElementsAlongPolyline const&)>
                get_cached_item = &BoundaryElementsAlongPolyline::getPolyline;
            return MeshGeoToolsLib::getBoundaryElements(
                _boundary_elements_along_polylines, get_cached_item,
                *dynamic_cast<const GeoLib::Polyline*>(&geoObj), _mesh,
                _mshNodeSearcher, false);
        }
        break;
        case GeoLib::GEOTYPE::SURFACE:
        {
            std::function<GeoLib::Surface(BoundaryElementsOnSurface const&)>
                get_cached_item = &BoundaryElementsOnSurface::getSurface;
            return MeshGeoToolsLib::getBoundaryElements(
                _boundary_elements_along_surfaces, get_cached_item,
                *dynamic_cast<const GeoLib::Surface*>(&geoObj), _mesh,
                _mshNodeSearcher, false);
        }
        break;
        default:
            const static std::vector<MeshLib::Element*> dummy(0);
            return dummy;
    }
}

}  // end namespace MeshGeoToolsLib
