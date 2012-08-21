/*
 * GMSHPolygonTree.h
 *
 *  Created on: Mar 27, 2012
 *      Author: fischeth
 */

#ifndef GMSHPOLYGONTREE_H_
#define GMSHPOLYGONTREE_H_

#include <vector>
#include <string>

// FileIO
#include "GMSHMeshDensityStrategy.h"
#include "GMSHPoint.h"
#include "GMSHLine.h"

// GEOLIB
#include "SimplePolygonTree.h"
#include "GEOObjects.h"
#include "PolylineWithSegmentMarker.h"

namespace FileIO {

class GMSHPolygonTree: public GEOLIB::SimplePolygonTree {
public:
	GMSHPolygonTree(GEOLIB::Polygon* polygon, GMSHPolygonTree * parent,
					GEOLIB::GEOObjects &geo_objs, std::string const& geo_name,
					GMSHMeshDensityStrategy * mesh_density_strategy);
	virtual ~GMSHPolygonTree();

	/**
	 * If the station point is inside the polygon, the method inserts the station into
	 * the internal vector of stations. This method works recursive!
	 * @param pnt the station point
	 * @return true if the station is inside the polygon
	 */
	bool insertStation(GEOLIB::Point const* pnt);
	/**
	 * If at least one (end) point (of a line segment) of the polyline is inside the polygon
	 * the polyline is inserted to the internal vector of polylines.
	 *
	 * Intersection points are inserted into the points vector the polygon and the polyline
	 * are based on. The id of the intersection point is inserted in both the polygon and the
	 * polyline, i.e. the two intersecting line segments are splitt into four line segment.
	 *
	 * Line segments of the polyline that are completely within the polygon are inserted into
	 * the internal vector _gmsh_lines_for_constraints. The childs of this GMSHPolygonTree node
	 * are checked recursively.
	 * @param ply the polyline that should be inserted
	 */
	void insertPolyline(GEOLIB::PolylineWithSegmentMarker * ply);

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
	void createGMSHPoints(std::vector<FileIO::GMSHPoint*> & gmsh_pnts) const;

	virtual void writeLineLoop(size_t &line_offset, size_t &sfc_offset, std::ostream& out) const;
	void writeSubPolygonsAsLineConstraints(size_t &line_offset, size_t sfc_number, std::ostream& out) const;
	virtual void writeLineConstraints(size_t &line_offset, size_t sfc_number, std::ostream& out) const;
	void writeStations(size_t & pnt_id_offset, size_t sfc_number, std::ostream& out) const;
	void writeAdditionalPointData(size_t & pnt_id_offset, size_t sfc_number, std::ostream& out) const;

private:
	void getPointsFromSubPolygons(std::vector<GEOLIB::Point const*>& pnts);
	void getStationsInsideSubPolygons(std::vector<GEOLIB::Point const*>& stations);
	const std::list<SimplePolygonTree*>& getChilds() const;
	const std::list<GEOLIB::GEOObjects*>& getGeoObjects () const;

	GEOLIB::GEOObjects & _geo_objs;
	std::string const& _geo_name;
	std::vector<GEOLIB::Point const*> _stations;
	std::vector<GEOLIB::PolylineWithSegmentMarker*> _plys;
	std::vector<FileIO::GMSHLine*> _gmsh_lines_for_constraints;

	GMSHMeshDensityStrategy * _mesh_density_strategy;
};

}

#endif /* GMSHPOLYGONTREE_H_ */
