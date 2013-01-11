/**
 * \file
 * \author Thomas Fischer
 * \date   Mar 27, 2012
 * \brief  Implementation of the GMSHPolygonTree class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshIO/GMSHPolygonTree.h"

#include "GMSHNoMeshDensity.h"
#include "GMSHFixedMeshDensity.h"
#include "GMSHAdaptiveMeshDensity.h"

namespace FileIO {

GMSHPolygonTree::GMSHPolygonTree(GeoLib::Polygon* polygon, GMSHPolygonTree* parent,
				GeoLib::GEOObjects &geo_objs, std::string const& geo_name,
				GMSHMeshDensityStrategy * mesh_density_strategy) :
	GeoLib::SimplePolygonTree(polygon, parent), _geo_objs(geo_objs), _geo_name(geo_name),
	_mesh_density_strategy(mesh_density_strategy)
{}

GMSHPolygonTree::~GMSHPolygonTree()
{}

bool GMSHPolygonTree::insertStation(GeoLib::Point const* station)
{

	if (_node_polygon->isPntInPolygon(*station)) {
		// try to insert station into the child nodes
		for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin());
			 it != _childs.end(); it++) {
			if (((*it)->getPolygon())->isPntInPolygon (*station)) {
				return dynamic_cast<GMSHPolygonTree*>((*it))->insertStation (station);
			}
		}
		// station did not fit into child nodes -> insert the station into this node
		_stations.push_back (station);
		return true;
	} else {
		return false;
	}
}

void GMSHPolygonTree::insertPolyline (GeoLib::PolylineWithSegmentMarker * ply)
{
	if (_node_polygon->isPartOfPolylineInPolygon(*ply)) {
		// check childs
		for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin());
			it != _childs.end(); it++) {
			dynamic_cast<GMSHPolygonTree*>((*it))->insertPolyline (ply);
		}
		_plys.push_back(ply);

		// calculate possible intersection points
		// pay attention: loop bound is not fix!
		size_t n_segments (ply->getNumberOfPoints()-1);
		GeoLib::Point tmp_pnt;
		for (size_t k(0); k<n_segments; k++) {
			if (! ply->isSegmentMarked(k)) {
				size_t seg_num(0);
				GeoLib::Point *intersection_pnt(new GeoLib::Point);
				while (_node_polygon->getNextIntersectionPointPolygonLine(*(ply->getPoint(k)), *(ply->getPoint(k+1)), intersection_pnt, seg_num)) {
					// insert the intersection point to point vector of GEOObjects instance
					size_t pnt_id;
					// appendPoints deletes an already existing point
					if (_geo_objs.appendPoint(intersection_pnt, _geo_name, pnt_id)) {
						// case: new point
						// modify the polygon
						_node_polygon->insertPoint(seg_num+1, pnt_id);
						// modify the polyline
						ply->insertPoint(k+1, pnt_id);
						n_segments++;
					} else {
						// case: existing point
						// check if point id is within the polygon
						if (! _node_polygon->isPointIDInPolyline(pnt_id)) {
							_node_polygon->insertPoint(seg_num+1, pnt_id);
						}
						// check if point id is in polyline
						if (! ply->isPointIDInPolyline(pnt_id)) {
							ply->insertPoint(k+1, pnt_id);
							n_segments++;
						}
					}

					size_t tmp_seg_num(seg_num+1);
					if (! _node_polygon->getNextIntersectionPointPolygonLine(*(ply->getPoint(k)), *(ply->getPoint(k+1)), &tmp_pnt, tmp_seg_num)) {
						// check a point of the segment except the end points
						for (size_t i(0); i<3; i++) {
							tmp_pnt[i] = ((*(ply->getPoint(k)))[i] + (*(ply->getPoint(k+1)))[i]) / 2;
						}
						if (_node_polygon->isPntInPolygon(tmp_pnt)) {
							ply->markSegment(k, true);
							// insert line segment as constraint
							_gmsh_lines_for_constraints.push_back(new GMSHLine(ply->getPointID(k), ply->getPointID(k+1)));
						}
					}
					intersection_pnt = new GeoLib::Point;
					seg_num++;
				}

				// check a point of the segment except the end points
				for (size_t i(0); i<3; i++) {
					tmp_pnt[i] = ((*(ply->getPoint(k)))[i] + (*(ply->getPoint(k+1)))[i]) / 2;
				}
				if (_node_polygon->isPntInPolygon(tmp_pnt)) {
					ply->markSegment(k, true);
					// insert line segment as constraint
					_gmsh_lines_for_constraints.push_back(new GMSHLine(ply->getPointID(k), ply->getPointID(k+1)));
				}
			}
		}
	}
}

void GMSHPolygonTree::initMeshDensityStrategy()
{
	if (dynamic_cast<GMSHAdaptiveMeshDensity*> (_mesh_density_strategy)) {
		// collect points
		std::vector<GeoLib::Point const*> pnts;
		const size_t n_pnts_polygon (_node_polygon->getNumberOfPoints());
		for (size_t k(0); k<n_pnts_polygon; k++) {
			pnts.push_back(_node_polygon->getPoint(k));
		}
		getPointsFromSubPolygons(pnts);

		const size_t n_plys (_plys.size());
		for (size_t k(0); k<n_plys; k++) {
			const size_t n_pnts_in_kth_ply(_plys[k]->getNumberOfPoints());
			for (size_t j(0); j<n_pnts_in_kth_ply; j++) {
				pnts.push_back(_plys[k]->getPoint(j));
			}
		}

		// give collected points to the mesh density strategy
		_mesh_density_strategy->init(pnts);
		// insert constraints
		dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->addPoints(_stations);
		std::vector<GeoLib::Point const*> stations;
		getStationsInsideSubPolygons(stations);
		dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->addPoints(stations);
	}
}

void GMSHPolygonTree::createGMSHPoints(std::vector<FileIO::GMSHPoint*> & gmsh_pnts) const
{
	const size_t n_pnts_polygon (_node_polygon->getNumberOfPoints());
	for (size_t k(0); k<n_pnts_polygon; k++) {
		const size_t id (_node_polygon->getPointID(k));
		GeoLib::Point const*const pnt(_node_polygon->getPoint(k));
		gmsh_pnts[id] = new GMSHPoint(*pnt, id, _mesh_density_strategy->getMeshDensityAtPoint(pnt));
	}

	const size_t n_plys(_plys.size());
	for (size_t k(0); k<n_plys; k++) {
		const size_t n_pnts_in_ply(_plys[k]->getNumberOfPoints());
		for (size_t j(0); j<n_pnts_in_ply; j++) {
			if (_node_polygon->isPntInPolygon(*(_plys[k]->getPoint(j)))) {
				const size_t id (_plys[k]->getPointID(j));
				GeoLib::Point const*const pnt(_plys[k]->getPoint(j));
				gmsh_pnts[id] = new GMSHPoint(*pnt, id, _mesh_density_strategy->getMeshDensityAtPoint(pnt));
			}
		}
	}

	// walk through childs
	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin()); it != _childs.end(); it++) {
		dynamic_cast<GMSHPolygonTree*>((*it))->createGMSHPoints(gmsh_pnts);
	}
}

void GMSHPolygonTree::writeLineLoop(size_t &line_offset, size_t &sfc_offset, std::ostream& out) const
{
	const size_t n_pnts (_node_polygon->getNumberOfPoints());
	size_t first_pnt_id(_node_polygon->getPointID(0)), second_pnt_id;
	for (size_t k(1); k<n_pnts; k++) {
		second_pnt_id = _node_polygon->getPointID(k);
		out << "Line(" << line_offset + k-1 << ") = {" << first_pnt_id << "," << second_pnt_id << "};" << std::endl;
		first_pnt_id = second_pnt_id;
	}
	out << "Line Loop(" << line_offset + n_pnts-1 << ") = {";
	for (size_t k(0); k<n_pnts - 2; k++) {
		out << line_offset+k << ",";
	}
	out << line_offset+n_pnts-2 << "};" << std::endl;
	out << "Plane Surface(" << sfc_offset << ") = {" << line_offset+n_pnts-1 << "};" << std::endl;
	line_offset += n_pnts;
	sfc_offset++;
}

void GMSHPolygonTree::writeLineConstraints(size_t &line_offset, size_t sfc_number, std::ostream& out) const
{
	const size_t n_plys (_plys.size());
	for (size_t j(0); j<n_plys; j++) {
		const size_t n_pnts(_plys[j]->getNumberOfPoints());
		size_t first_pnt_id(_plys[j]->getPointID(0)), second_pnt_id;
		for (size_t k(1); k<n_pnts; k++) {
			second_pnt_id = _plys[j]->getPointID(k);
			if (_plys[j]->isSegmentMarked(k-1) && _node_polygon->isPntInPolygon(*(_plys[j]->getPoint(k)))) {
				out << "Line(" << line_offset + k-1 << ") = {" << first_pnt_id << "," << second_pnt_id << "};" << std::endl;
				out << "Line { " << line_offset+k-1 << " } In Surface { " << sfc_number << " };" << std::endl;
			}
			first_pnt_id = second_pnt_id;
		}
		line_offset += n_pnts;
	}
}

void GMSHPolygonTree::writeSubPolygonsAsLineConstraints(size_t &line_offset, size_t sfc_number, std::ostream& out) const
{
	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin()); it != _childs.end(); it++) {
		dynamic_cast<GMSHPolygonTree*>((*it))->writeSubPolygonsAsLineConstraints(line_offset, sfc_number, out);
	}

	if (_parent != NULL) {
		const size_t n_pnts(_node_polygon->getNumberOfPoints());
		size_t first_pnt_id(_node_polygon->getPointID(0)), second_pnt_id;
		for (size_t k(1); k<n_pnts; k++) {
			second_pnt_id = _node_polygon->getPointID(k);
			out << "Line(" << line_offset + k-1 << ") = {" << first_pnt_id << "," << second_pnt_id << "};" << std::endl;
			first_pnt_id = second_pnt_id;
			out << "Line { " << line_offset+k-1 << " } In Surface { " << sfc_number << " };" << std::endl;
		}
		line_offset += n_pnts;
	}

}

void GMSHPolygonTree::writeStations(size_t & pnt_id_offset, size_t sfc_number, std::ostream& out) const
{
	const size_t n_stations(_stations.size());
	for (size_t k(0); k<n_stations; k++) {
		out << "Point(" << pnt_id_offset + k << ") = {" << (*(_stations[k]))[0] << "," << (*(_stations[k]))[1] << ", 0.0, ";
		out << _mesh_density_strategy->getMeshDensityAtPoint(_stations[k]) << "};" << std::endl;
		out << "Point { " << pnt_id_offset + k << " } In Surface { " << sfc_number << " }; " << std::endl;
	}
	pnt_id_offset += n_stations;
}

void GMSHPolygonTree::writeAdditionalPointData(size_t & pnt_id_offset, size_t sfc_number, std::ostream& out) const
{
	if (dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)) {
		std::vector<GeoLib::Point*> steiner_pnts;
		dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->getSteinerPoints(steiner_pnts, 0);
		const size_t n(steiner_pnts.size());
		for (size_t k(0); k<n; k++) {
			if (_node_polygon->isPntInPolygon(*(steiner_pnts[k]))) {
				out << "Point(" << pnt_id_offset + k << ") = {" << (*(steiner_pnts[k]))[0] << "," << (*(steiner_pnts[k]))[1] << ", 0.0, ";
				out << _mesh_density_strategy->getMeshDensityAtPoint(steiner_pnts[k]) << "};" << std::endl;
				out << "Point { " << pnt_id_offset + k << " } In Surface { " << sfc_number << " }; " << std::endl;
			}
			delete steiner_pnts[k];
		}
		pnt_id_offset += n;
	}

#ifndef NDEBUG
	if (dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)) {
		std::vector<GeoLib::Point*> * pnts(new std::vector<GeoLib::Point*>);
		std::vector<GeoLib::Polyline*> *plys(new std::vector<GeoLib::Polyline*>);
		dynamic_cast<GMSHAdaptiveMeshDensity*>(_mesh_density_strategy)->getQuadTreeGeometry(*pnts, *plys);
		std::string quad_tree_geo("QuadTree");
		_geo_objs.addPointVec(pnts, quad_tree_geo);
		std::vector<size_t> const& id_map ((_geo_objs.getPointVecObj(quad_tree_geo))->getIDMap());
		for (size_t k(0); k<plys->size(); k++) {
			for (size_t j(0); j<(*plys)[k]->getNumberOfPoints(); j++) {
				((*plys)[k])->setPointID(j, id_map[((*plys)[k])->getPointID(j)]);
			}
		}
		_geo_objs.addPolylineVec(plys, quad_tree_geo);
	}
#endif

}

void GMSHPolygonTree::getPointsFromSubPolygons(std::vector<GeoLib::Point const*>& pnts)
{
	const size_t n_pnts_polygon (_node_polygon->getNumberOfPoints());
	for (size_t k(0); k<n_pnts_polygon; k++) {
		pnts.push_back(_node_polygon->getPoint(k));
	}

	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin()); it != _childs.end(); it++) {
		dynamic_cast<GMSHPolygonTree*>((*it))->getPointsFromSubPolygons(pnts);
	}
}

void GMSHPolygonTree::getStationsInsideSubPolygons(std::vector<GeoLib::Point const*>& stations)
{
	const size_t n_stations(_stations.size());
	for (size_t k(0); k<n_stations; k++) {
		stations.push_back(_stations[k]);
	}

	for (std::list<SimplePolygonTree*>::const_iterator it (_childs.begin()); it != _childs.end(); it++) {
		dynamic_cast<GMSHPolygonTree*>((*it))->getStationsInsideSubPolygons(stations);
	}
}

} // end namespace FileIO
