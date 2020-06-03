/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <vector>

namespace GeoLib
{
struct GeoObject;
class Point;
class Polyline;
class Surface;
}

namespace MeshLib
{
class Mesh;
class Element;
}

namespace MeshGeoToolsLib
{
class MeshNodeSearcher;
class BoundaryElementsAtPoint;
class BoundaryElementsAlongPolyline;
class BoundaryElementsOnSurface;

/**
 * This class searches boundary elements located on a given geometric object, i.e. polyline and surface.
 * Note that internal boundaries are currently not supported.
 */
class BoundaryElementsSearcher
{
public:
    /**
     * Constructor
     * @param mesh             a mesh object
     * @param mshNodeSearcher  a MeshNodeSearcher object which is internally used to search mesh nodes
     */
    BoundaryElementsSearcher(MeshLib::Mesh const& mesh,
                             MeshNodeSearcher const& mshNodeSearcher);

    /// destructor
    virtual ~BoundaryElementsSearcher();

    /**
     * generate boundary elements on the given geometric object (point, polyline, surface).
     *
     * @param geoObj a GeoLib::GeoObject where the nearest mesh node is searched for
     * @param multiple_nodes_allowed allows for finding multiple nodes within
     * the given search radius for a point
     * @return a vector of boundary element objects
     */
    std::vector<MeshLib::Element*> const& getBoundaryElements(
        GeoLib::GeoObject const& geoObj, bool const multiple_nodes_allowed);

    /**
     * generate boundary elements at the given point.
     * @param point Search the mesh for given point
     * @param multiple_nodes_allowed allows for finding multiple nodes within
     * the given search radius
     * @return a vector of boundary elements
     */
    std::vector<MeshLib::Element*> const& getBoundaryElementsAtPoint(
        GeoLib::Point const& point, bool const multiple_nodes_allowed);

    /**
     * generate boundary elements on the given polyline.
     * @param polyline the GeoLib::Polyline the nearest mesh nodes are searched
     * for
     * @return a vector of boundary element objects
     */
    std::vector<MeshLib::Element*> const& getBoundaryElementsAlongPolyline(
        GeoLib::Polyline const& polyline);

    /**
     * generate boundary elements on the given surface.
     * @param surface the GeoLib::Surface the nearest mesh nodes are searched
     * for
     * @return a vector of boundary element objects
     */
    std::vector<MeshLib::Element*> const& getBoundaryElementsOnSurface(
        GeoLib::Surface const& surface);

private:
    MeshLib::Mesh const& mesh_;
    MeshNodeSearcher const& mshNodeSearcher_;
    std::vector<BoundaryElementsAtPoint*> boundary_elements_at_point_;
    std::vector<BoundaryElementsAlongPolyline*> boundary_elements_along_polylines_;
    std::vector<BoundaryElementsOnSurface*> boundary_elements_along_surfaces_;
};

} // end namespace MeshGeoToolsLib
