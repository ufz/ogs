/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-25
 * \brief  Implementation of the GeoMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "GeoMapper.h"

#include <numeric>

#include "FileIO/AsciiRasterInterface.h"

#include "AABB.h"
#include "Mesh.h"
#include "Elements/Element.h"
#include "Node.h"
#include "MeshSurfaceExtraction.h"
#include "AnalyticalGeometry.h"
#include "PointWithID.h"
#include "Raster.h"
#include "readMeshFromFile.h"
#include "StationBorehole.h"


GeoMapper::GeoMapper(GeoLib::GEOObjects &geo_objects, const std::string &geo_name)
	: _geo_objects(geo_objects), _geo_name(const_cast<std::string&>(geo_name)), _mesh(nullptr), _grid(nullptr), _raster(nullptr)
{
}

GeoMapper::~GeoMapper()
{
	delete _mesh;
	delete _raster;
}

void GeoMapper::mapOnDEM(const std::string &file_name)
{
	this->_raster = FileIO::AsciiRasterInterface::getRasterFromASCFile(file_name);
	if (! _raster) {
		ERR("GeoMapper::mapOnDEM(): failed to load %s", file_name.c_str());
		return;
	}
	this->mapData();
}

void GeoMapper::mapOnMesh(const std::string &file_name)
{
	MeshLib::Mesh *mesh (FileIO::readMeshFromFile(file_name));
	mapOnMesh(mesh);
	delete mesh;
}

void GeoMapper::mapOnMesh(const MeshLib::Mesh* mesh)
{
	if (mesh->getDimension()<3)
		this->_mesh = new MeshLib::Mesh(*mesh);
	else
	{
		const MathLib::Vector3 dir(0,0,-1);
		this->_mesh = MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir);
	}
	std::vector<GeoLib::PointWithID*> sfc_pnts;
	// init grid
	_grid = this->getFlatGrid(_mesh, sfc_pnts);
	this->mapData();

	delete _grid;

	const size_t n_sfc_pnts(sfc_pnts.size());
	for (size_t k(0); k<n_sfc_pnts; k++) {
		delete sfc_pnts[k];
	}
}

void GeoMapper::mapData()
{
	const std::vector<GeoLib::Point*> *points (this->_geo_objects.getPointVec(this->_geo_name));
	GeoLib::Station* stn_test = dynamic_cast<GeoLib::Station*>((*points)[0]);
	bool is_borehole(false);
	if (stn_test != nullptr && static_cast<GeoLib::StationBorehole*>((*points)[0])->type() == GeoLib::Station::StationType::BOREHOLE)
		is_borehole = true;
	
	double min_val(0), max_val(0);
	if (_mesh)
	{
		GeoLib::AABB<GeoLib::Point> bounding_box (_mesh->getNodes().begin(), _mesh->getNodes().end());
		min_val = bounding_box.getMinPoint()[2];
		max_val = bounding_box.getMaxPoint()[2];
	}
	size_t nPoints (points->size());

	if (!is_borehole)
	{
		for (unsigned j=0; j<nPoints; ++j)
		{
			GeoLib::Point* pnt ((*points)[j]);
			(*pnt)[2] = (_grid) ? this->getMeshElevation((*pnt)[0],(*pnt)[1], min_val, max_val)
				                : this->getDemElevation(*pnt);
		}
	}
	else
	{
		for (unsigned j=0; j<nPoints; ++j)
		{
			GeoLib::Point* pnt ((*points)[j]);
			double offset = (_grid) ? (this->getMeshElevation((*pnt)[0],(*pnt)[1], min_val, max_val) - (*pnt)[2])
				                    :  this->getDemElevation(*pnt);

			GeoLib::StationBorehole* borehole = static_cast<GeoLib::StationBorehole*>(pnt);
			const std::vector<GeoLib::Point*> layers = borehole->getProfile();
			size_t nLayers = layers.size();
			for (unsigned k=0; k<nLayers; ++k)
			{
				GeoLib::Point* layer_pnt = layers[k];
				(*layer_pnt)[2] = (*layer_pnt)[2] + offset;
			}
		}
	}
}

float GeoMapper::getDemElevation(GeoLib::Point const& pnt) const
{
	double const elevation (_raster->getValueAtPoint(pnt));
	if (std::abs(elevation-_raster->getNoDataValue()) < std::numeric_limits<double>::epsilon())
		return 0.0;
	return static_cast<float>(elevation);
}

double GeoMapper::getMeshElevation(double x, double y, double min_val, double max_val) const
{
	double coords[3] = {x,y,0};
	const GeoLib::PointWithID* pnt = _grid->getNearestPoint(coords);
	const std::vector<MeshLib::Element*> elements (_mesh->getNode(pnt->getID())->getElements());
	GeoLib::Point* intersection (nullptr);

	for (std::size_t i=0; i<elements.size(); ++i)
	{
		if (intersection==nullptr && elements[i]->getGeomType() != MeshElemType::LINE)
			intersection=GeoLib::triangleLineIntersection(*elements[i]->getNode(0), *elements[i]->getNode(1), *elements[i]->getNode(2), GeoLib::Point(x,y,max_val), GeoLib::Point(x,y,min_val));
		if (intersection==nullptr && elements[i]->getGeomType() == MeshElemType::QUAD)
			intersection=GeoLib::triangleLineIntersection(*elements[i]->getNode(0), *elements[i]->getNode(2), *elements[i]->getNode(3), GeoLib::Point(x,y,max_val), GeoLib::Point(x,y,min_val));
	}
	if (intersection)
		return (*intersection)[2];
	// if something goes wrong, simply take the elevation of the nearest mesh node
	return (*(_mesh->getNode(pnt->getID())))[2];
}

GeoLib::Grid<GeoLib::PointWithID>* GeoMapper::getFlatGrid(MeshLib::Mesh const*const mesh, std::vector<GeoLib::PointWithID*> sfc_pnts) const
{
	size_t nNodes (mesh->getNNodes());
	sfc_pnts.resize(nNodes);
	const std::vector<MeshLib::Node*> nodes (mesh->getNodes());
	for (unsigned i(0); i<nNodes; ++i)
		sfc_pnts[i] = new GeoLib::PointWithID((*nodes[i])[0], (*nodes[i])[1], 0.0, nodes[i]->getID());

	return new GeoLib::Grid<GeoLib::PointWithID>(sfc_pnts.begin(), sfc_pnts.end());
}

unsigned getIndexInPntVec(GeoLib::Point const*const pnt, std::vector<GeoLib::Point*> const*const points)
{
	auto it (std::find(points->begin(), points->end(), pnt));
	return static_cast<unsigned>(std::distance(points->begin(), it));
}

std::vector<GeoLib::Polyline*>* copyPolylinesVector(const std::vector<GeoLib::Polyline*> *polylines, std::vector<GeoLib::Point*> *points)
{
	std::size_t nLines = polylines->size();
	std::vector<GeoLib::Polyline*> *new_lines = new std::vector<GeoLib::Polyline*>(nLines);
	for (std::size_t i=0; i<nLines; ++i)
	{
		(*new_lines)[i] = new GeoLib::Polyline(*points);
		std::size_t nLinePnts ((*polylines)[i]->getNumberOfPoints());
		for (std::size_t j=0; j<nLinePnts; ++j)
			(*new_lines)[i]->addPoint((*polylines)[i]->getPointID(j));
	}
	return new_lines;
}


void GeoMapper::advancedMapOnMesh(const MeshLib::Mesh* mesh, const std::string &new_geo_name)
{
	const std::vector<GeoLib::Point*> *points (this->_geo_objects.getPointVec(this->_geo_name));
	const std::vector<GeoLib::Polyline*> *org_lines (this->_geo_objects.getPolylineVec(this->_geo_name));

	const GeoLib::AABB<GeoLib::Point> aabb(points->begin(), points->end());
	const double eps = sqrt(std::numeric_limits<float>::epsilon()) *
		               sqrt( MathLib::sqrDist(aabb.getMinPoint(),aabb.getMaxPoint())) ;

	// copy geometry (and set z=0 for all points)
	unsigned nGeoPoints ( points->size() );
	std::vector<GeoLib::Point*> *new_points = new std::vector<GeoLib::Point*>(nGeoPoints);
	for (size_t i=0; i<nGeoPoints; ++i)
		(*new_points)[i] = new GeoLib::Point((*(*points)[i])[0],(*(*points)[i])[1],0.0);
	std::vector<GeoLib::Polyline*> *new_lines (copyPolylinesVector(this->_geo_objects.getPolylineVec(this->_geo_name), new_points));

	GeoLib::Grid<GeoLib::Point> grid(new_points->begin(), new_points->end());
	double max_segment_length (this->getMaxSegmentLength(*new_lines));
	max_segment_length *= max_segment_length; // squared so it can be compared to the squared distances calculated later
	
	const unsigned nMeshNodes ( mesh->getNNodes() );	
	std::vector<int> closest_geo_point(nMeshNodes); // index of closest geo point for each mesh node in (x,y)-plane
	std::vector<double> dist(nMeshNodes);  // distance between geo points and mesh nodes in (x,y)-plane
	for (size_t i=0; i<nMeshNodes; ++i)
	{
		const double zero_coords[3] = {(* mesh->getNode(i))[0], (* mesh->getNode(i))[1], 0.0};
		GeoLib::Point* pnt = grid.getNearestPoint(zero_coords);
		dist[i] = MathLib::sqrDist(pnt->getCoords(), zero_coords);
		closest_geo_point[i] = (dist[i]<=max_segment_length) ? getIndexInPntVec(pnt, new_points) : -1;
	}
	
	// store for each point the line segment to which it was added.
	const size_t nLines (new_lines->size());
	std::vector< std::vector<unsigned> > line_segment_map(nLines);
	for (std::size_t i=0; i<nLines; ++i)
	{
		line_segment_map[i] = std::vector<unsigned>((*new_lines)[i]->getNumberOfPoints(),0);
		std::iota(line_segment_map[i].begin(), line_segment_map[i].end(), 0);
	}

	for (std::size_t i=0; i<nMeshNodes; ++i)
	{
		// if mesh node too far away or exactly at point position
		if (closest_geo_point[i] == -1 || dist[i] < eps) continue; 

		const MeshLib::Node* node (mesh->getNode(i));
		for (std::size_t l=0; l<nLines; ++l)
		{
			// find relevant polylines
			if (!(*org_lines)[l]->isPointIDInPolyline(closest_geo_point[i])) continue;
			
			// find point position of closest geo point in original polyline
			GeoLib::Polyline* ply ((*org_lines)[l]);
			std::size_t nLinePnts ( ply->getNumberOfPoints() );
			std::size_t node_index_in_ply (0);
			for (node_index_in_ply=0; node_index_in_ply<nLinePnts; ++node_index_in_ply)
				if (ply->getPoint(node_index_in_ply) == (*points)[closest_geo_point[i]])
					break;
			const GeoLib::Point* geo_point (ply->getPoint(node_index_in_ply));

			// check if line segments connected to closest geo point intersect connected elements of current node
			const std::vector<MeshLib::Element*> elements (node->getElements());
			const std::size_t nElems = elements.size();
			for (std::size_t e=0; e<nElems; ++e)
			{
				const unsigned nEdges (elements[e]->getNEdges());
				unsigned intersection_count (0);

				for (unsigned n=0; n<nEdges; ++n)
				{
					if (intersection_count>1) break; //already two intersections

					const MeshLib::Element* line = elements[e]->getEdge(n);
					unsigned index_offset(0); // default: add to first line segment
					GeoLib::Point* intersection (NULL);
					if (node_index_in_ply>0) // test line segment before closest point
						intersection = calcIntersection(line->getNode(0), line->getNode(1), geo_point, ply->getPoint(node_index_in_ply-1));
					if (intersection == NULL && node_index_in_ply<(nLinePnts-1)) // test line segment after closest point
					{
						intersection = calcIntersection(line->getNode(0), line->getNode(1), geo_point, ply->getPoint(node_index_in_ply+1));
						index_offset = 1; // add to second segment
					}
					if (intersection)
					{
						intersection_count++;
						unsigned start_point_idx = static_cast<unsigned>(std::distance(line_segment_map[l].begin(), std::find_if(line_segment_map[l].begin(), line_segment_map[l].end(), [&node_index_in_ply, &index_offset](unsigned a){return a==node_index_in_ply+index_offset-1;})));
						unsigned end_point_idx   = static_cast<unsigned>(std::distance(line_segment_map[l].begin(), std::find_if(line_segment_map[l].begin(), line_segment_map[l].end(), [&node_index_in_ply, &index_offset](unsigned a){return a==node_index_in_ply+index_offset;})));
						std::size_t pos = getPointPosInLine((*new_lines)[l], start_point_idx, end_point_idx, intersection, eps);

						if (pos)
						{
							const std::size_t pnt_pos (new_points->size());
							new_points->push_back(intersection);
							(*new_lines)[l]->insertPoint(pos, pnt_pos);
							line_segment_map[l].insert(line_segment_map[l].begin()+pos, node_index_in_ply+index_offset-1);
						}
					}
				}
			}
		}
	}

	this->_geo_objects.addPointVec(new_points, const_cast<std::string&>(new_geo_name));
	std::vector<size_t> pnt_id_map = this->_geo_objects.getPointVecObj(new_geo_name)->getIDMap();
	for (std::size_t i=0; i<new_lines->size(); ++i)
		(*new_lines)[i]->updatePointIDs(pnt_id_map);
	this->_geo_objects.addPolylineVec(new_lines, new_geo_name);

	// map new geometry incl. additional point using the normal mapping method
	this->_geo_name = new_geo_name;
	this->mapOnMesh(mesh);
}

GeoLib::Point* GeoMapper::calcIntersection(GeoLib::Point const*const p1, GeoLib::Point const*const p2, GeoLib::Point const*const q1, GeoLib::Point const*const q2) const
{
	const double x1 = (*p1)[0], x2 = (*p2)[0], x3 = (*q1)[0], x4 = (*q2)[0];
	const double y1 = (*p1)[1], y2 = (*p2)[1], y3 = (*q1)[1], y4 = (*q2)[1];
 
	const double det = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
	if (fabs(det) < std::numeric_limits<double>::epsilon()) return NULL;
 
	const double pre  = (x1*y2 - y1*x2);
	const double post = (x3*y4 - y3*x4);
	const double x = ( pre * (x3 - x4) - (x1 - x2) * post ) / det;
	const double y = ( pre * (y3 - y4) - (y1 - y2) * post ) / det;
 
	// Check if the x and y coordinates are within both line segments
	if (isPntInBoundingBox(x1,y1,x2,y2,x,y) && isPntInBoundingBox(x3,y3,x4,y4,x,y))
		return new GeoLib::Point(x, y, 0);
	return NULL;
}

unsigned GeoMapper::getPointPosInLine(GeoLib::Polyline const*const line, unsigned start, unsigned end, GeoLib::Point const*const point, double eps) const
{
	const double* first = line->getPoint(start)->getCoords();
	const double* pnt   = point->getCoords();
	const double max_dist = MathLib::sqrDist(first, pnt);

	// if point is at start or end of line segment
	if (max_dist<eps && MathLib::sqrDist(pnt, line->getPoint(end)->getCoords())) return 0;

	for (std::size_t i=start+1; i<end; ++i)
	{
		const double* current = (*line->getPoint(i)).getCoords();
		if (MathLib::sqrDist(pnt, current) < eps) return 0;
		if (MathLib::sqrDist(first, current) > max_dist) return i;
	}
	return end; // last point of segment
}

bool GeoMapper::isPntInBoundingBox(double ax, double ay, double bx, double by, double px, double py) const
{
	if ( px < (std::min(ax, bx)-std::numeric_limits<double>::epsilon()) || px > (std::max(ax, bx)+std::numeric_limits<double>::epsilon()) || 
		 py < (std::min(ay, by)-std::numeric_limits<double>::epsilon()) || py > (std::max(ay, by)+std::numeric_limits<double>::epsilon()) ) 
		 return false;
	return true;
}

double GeoMapper::getMaxSegmentLength(const std::vector<GeoLib::Polyline*> &lines) const
{
	double max_segment_length (0);
	const std::size_t nPlys ( lines.size() );
	for (size_t i=0; i<nPlys; ++i)
	{
		const GeoLib::Polyline* line = lines[i];
		const std::size_t nPlyPoints = line->getNumberOfPoints();
		for (size_t j=1; j<nPlyPoints; ++j)
		{
			const double dist (line->getLength(j)-line->getLength(j-1));
			if (dist>max_segment_length)
				max_segment_length=dist;
		}	
	}
	return max_segment_length;
}






