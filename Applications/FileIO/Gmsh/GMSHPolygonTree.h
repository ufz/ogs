/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include <string>

#include "GMSHMeshDensityStrategy.h"
#include "GMSHPoint.h"
#include "GMSHLine.h"

#include "GeoLib/SimplePolygonTree.h"

namespace GeoLib
{
class GEOObjects;
class Polygon;
class PolylineWithSegmentMarker;
class PolygonWithSegmentMarker;
}

namespace FileIO
{
namespace GMSH
{

class GMSHPolygonTree: public GeoLib::SimplePolygonTree {
public:
    GMSHPolygonTree(GeoLib::PolygonWithSegmentMarker* polygon, GMSHPolygonTree * parent,
                    GeoLib::GEOObjects &geo_objs, std::string const& geo_name,
                    GMSHMeshDensityStrategy * mesh_density_strategy);
    virtual ~GMSHPolygonTree();

    /** Mark the segments shared by several polygons. */
    void markSharedSegments();

    /**
     * If the station point is inside the polygon, the method inserts the station into
     * the internal vector of stations. This method works recursive!
     * @param pnt the station point
     * @return true if the station is inside the polygon
     */
    bool insertStation(GeoLib::Point const* pnt);
    /**
     * If at least one (end) point (of a line segment) of the polyline is inside the polygon
     * the polyline is inserted to the internal vector of polylines.
     *
     * Intersection points are inserted into the points vector the polygon and the polyline
     * are based on. The id of the intersection point is inserted in both the polygon and the
     * polyline, i.e. the two intersecting line segments are splitt into four line segment.
     *
     * Line segments of the polyline that are completely within the polygon are inserted into
     * the internal vector _gmsh_lines_for_constraints. The children of this GMSHPolygonTree node
     * are checked recursively.
     * @param ply the polyline that should be inserted
     */
    void insertPolyline(GeoLib::PolylineWithSegmentMarker * ply);

    /**
     * Initialize the mesh density strategy with data. In case of GMSHAdaptiveMeshDensity
     * an instance of class QuadTree will be set up with the points within the top
     * level bounding polygon.
     */
    void initMeshDensityStrategy();

    /**
     * Method creates the gmsh point data structures - including the mesh density.
     * @param gmsh_pnts a vector of pointers to instances of class GMSHPoint
     */
    void createGMSHPoints(std::vector<GMSHPoint*> & gmsh_pnts) const;

    virtual void writeLineLoop(std::size_t &line_offset, std::size_t &sfc_offset, std::ostream& out) const;
    void writeSubPolygonsAsLineConstraints(std::size_t &line_offset, std::size_t sfc_number, std::ostream& out) const;
    virtual void writeLineConstraints(std::size_t &line_offset, std::size_t sfc_number, std::ostream& out) const;
    void writeStations(std::size_t & pnt_id_offset, std::size_t sfc_number, std::ostream& out) const;
    void writeAdditionalPointData(std::size_t & pnt_id_offset, std::size_t sfc_number, std::ostream& out) const;

private:
    void getPointsFromSubPolygons(std::vector<GeoLib::Point const*>& pnts);
    void getStationsInsideSubPolygons(std::vector<GeoLib::Point const*>& stations);
    void checkIntersectionsSegmentExistingPolylines(
        GeoLib::PolylineWithSegmentMarker* ply,
        GeoLib::Polyline::SegmentIterator const& segment_iterator);

    GeoLib::GEOObjects & _geo_objs;
    std::string const& _geo_name;
    std::vector<GeoLib::Point const*> _stations;
    std::vector<GeoLib::PolylineWithSegmentMarker*> _plys;
    std::vector<GMSHLine*> _gmsh_lines_for_constraints;

    GMSHMeshDensityStrategy * _mesh_density_strategy;
};

}  // end namespace GMSH
}  // end namespace FileIO
