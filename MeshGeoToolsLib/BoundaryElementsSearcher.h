/**
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */
#pragma once

#include <vector>

namespace GeoLib
{
class GeoObject;
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
     * @return a vector of boundary element objects
     */
    std::vector<MeshLib::Element*> const& getBoundaryElements(GeoLib::GeoObject const& geoObj);

    /// For each boundary element found the corresponding bulk element id and
    /// bulk element face id are returned.
    ///
    /// @param geometry The boundary elements, bulk element ids and  bulk
    /// element face ids will be extracted for the given geometry (point,
    /// polyline or surface).
    /// @return @copydoc BoundaryElementsAtPoint::_bulk_ids
    std::vector<std::pair<std::size_t, unsigned>> const& getBulkIDs(
        GeoLib::GeoObject const& geometry);

    /**
     * generate boundary elements at the given point.
     * @param point Search the mesh for given point
     * @return a vector of boundary elements
     */
    std::vector<MeshLib::Element*> const& getBoundaryElementsAtPoint(
        GeoLib::Point const& point);
    /// @copybrief BoundaryElementsSearcher::getBulkIDs()
    /// @param point The given point the boundary element will be searched for.
    /// @return @copydoc BoundaryElementsAtPoint::_bulk_ids
    std::vector<std::pair<std::size_t, unsigned>> const& getBulkIDsAtPoint(
        GeoLib::Point const& point);

    /**
     * generate boundary elements on the given polyline.
     * @param ply the GeoLib::Polyline the nearest mesh nodes are searched for
     * @return a vector of boundary element objects
     */
    std::vector<MeshLib::Element*> const& getBoundaryElementsAlongPolyline(
        GeoLib::Polyline const& ply);
    /// @copybrief BoundaryElementsSearcher::getBulkIDs()
    /// @param ply The given polyline the boundary element will be searched for.
    /// @return @copydoc BoundaryElementsAlongPolyline::_bulk_ids
    std::vector<std::pair<std::size_t, unsigned>> const&
    getBulkIDsAlongPolyline(GeoLib::Polyline const& ply);

    /**
     * generate boundary elements on the given surface.
     * @param sfc the GeoLib::Surface the nearest mesh nodes are searched for
     * @return a vector of boundary element objects
     */
    std::vector<MeshLib::Element*> const& getBoundaryElementsOnSurface(
        GeoLib::Surface const& sfc);
    /// @copybrief BoundaryElementsSearcher::getBulkIDs()
    /// @param sfc The given surface the boundary element will be searched
    /// for.
    /// @return @copydoc BoundaryElementsOnSurface::_bulk_ids
    std::vector<std::pair<std::size_t, unsigned>> const&
    getBulkIDsOnSurface(GeoLib::Surface const& sfc);

private:
    MeshLib::Mesh const& _mesh;
    MeshNodeSearcher const& _mshNodeSearcher;
    std::vector<BoundaryElementsAtPoint*> _boundary_elements_at_point;
    std::vector<BoundaryElementsAlongPolyline*> _boundary_elements_along_polylines;
    std::vector<BoundaryElementsOnSurface*> _boundary_elements_along_surfaces;
};

} // end namespace MeshGeoToolsLib
