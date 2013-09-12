/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-25
 * \brief  Implementation of the GeoMapper class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "GeoMapper.h"

#include <numeric>

#include "Mesh.h"
#include "Elements\Element.h"
#include "Elements\Tri.h"
#include "Node.h"
#include "MeshSurfaceExtraction.h"
#include "PointWithID.h"
#include "Raster.h"
#include "readMeshFromFile.h"
#include "StationBorehole.h"


GeoMapper::GeoMapper(GeoLib::GEOObjects &geo_objects, const std::string &geo_name)
	: _geo_objects(geo_objects), _geo_name(geo_name), _grid(NULL), _raster(nullptr)
{
}

GeoMapper::~GeoMapper()
{
	delete _raster;
}

void GeoMapper::mapOnDEM(const std::string &file_name)
{
	this->_raster = GeoLib::Raster::getRasterFromASCFile(file_name);
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
	std::vector<GeoLib::PointWithID*> sfc_pnts;
	// init grid
	_grid = this->getFlatGrid(mesh, sfc_pnts);
	this->mapData(mesh);

	delete _grid;

	const size_t n_sfc_pnts(sfc_pnts.size());
	for (size_t k(0); k<n_sfc_pnts; k++) {
		delete sfc_pnts[k];
	}
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


void GeoMapper::advancedMapOnMesh(const MeshLib::Mesh* mesh, std::string &new_geo_name)
{
	const double eps = sqrt(sqrt(std::numeric_limits<float>::epsilon()));

	// copy geometry (and set z=0 for all points)
	const std::vector<GeoLib::Point*> *points (this->_geo_objects.getPointVec(this->_geo_name));
	const std::vector<GeoLib::Polyline*> *org_lines (this->_geo_objects.getPolylineVec(this->_geo_name));
	const unsigned nGeoPoints ( points->size() );
	std::vector<GeoLib::Point*> *new_points = new std::vector<GeoLib::Point*>(nGeoPoints);
	for (size_t i=0; i<nGeoPoints; ++i)
		(*new_points)[i] = new GeoLib::Point((*(*points)[i])[0],(*(*points)[i])[1],0.0);
	std::vector<GeoLib::Polyline*> *new_lines (copyPolylinesVector(this->_geo_objects.getPolylineVec(this->_geo_name), new_points));

	GeoLib::Grid<GeoLib::Point> grid(new_points->begin(), new_points->end());
	double max_segment_length (this->getMaxSegmentLength(*new_lines));
	max_segment_length *= max_segment_length; // squared so it can be compared to the squared distances calculated later
	
	const unsigned nMeshNodes ( mesh->getNNodes() );	
	std::vector<int> closest_geo_point(nMeshNodes); 
	std::vector<double> dist(nMeshNodes); 
	for (size_t i=0; i<nMeshNodes; ++i)
	{
		const double* node_coords = mesh->getNode(i)->getCoords();
		const double zero_coords[3] = {node_coords[0], node_coords[1], 0.0};
		GeoLib::Point* pnt = grid.getNearestPoint(zero_coords);
		dist[i] = MathLib::sqrDist(pnt->getCoords(), zero_coords);
		if (dist[i]<=max_segment_length)
		{
			for (size_t j=0; j<nGeoPoints; ++j)
				if (pnt == (*new_points)[j])
					closest_geo_point[i] = j;
		}
		else
			closest_geo_point[i] = -1;
	}
	
	const size_t nLines (new_lines->size());
	// stores for each new point the line segment to which it was added.
	std::vector< std::vector<unsigned> > line_segment_map(nLines);
	for (std::size_t i=0; i<nLines; ++i)
	{
		line_segment_map[i] = std::vector<unsigned>((*new_lines)[i]->getNumberOfPoints(),0);
		std::iota(line_segment_map[i].begin(), line_segment_map[i].end(), 0);
	}

	for (std::size_t i=0; i<nMeshNodes; ++i)
	{
		if (closest_geo_point[i] == -1) continue; // is mesh node theoretically close enough?
		const MeshLib::Node* node (mesh->getNode(i));
		if (dist[i] < eps) // is mesh node == geo point?
		{
			(*(*new_points)[closest_geo_point[i]])[2] = (*node)[2];
			continue;
		}

		for (std::size_t l=0; l<nLines; ++l)
		{
			// find relevant polylines
			if (!(*org_lines)[l]->isPointIDInPolyline(closest_geo_point[i])) continue;
			
			GeoLib::Polyline* ply ((*org_lines)[l]);
			std::size_t nLinePnts ( ply->getNumberOfPoints() );
			std::size_t node_index_in_ply (0);
			for (node_index_in_ply=0; node_index_in_ply<nLinePnts; ++node_index_in_ply)
				if (ply->getPoint(node_index_in_ply) == (*points)[closest_geo_point[i]])
					break;

			const GeoLib::Point* geo_point (ply->getPoint(node_index_in_ply));

			// check if line intersects connected elements of current node
			const std::vector<MeshLib::Element*> elements (node->getElements());
			const std::size_t nElems = elements.size();
			for (std::size_t e=0; e<nElems; ++e)
			{
				const unsigned nNodes (elements[e]->getNNodes());
				const unsigned nEdges (elements[e]->getNEdges());
				unsigned intersection_count (0);

				for (unsigned n=0; n<nEdges; ++n)
				{
					if (intersection_count>1) break; //already two intersections

					const MeshLib::Element* line = elements[e]->getEdge(n);
					unsigned index_offset(0); // default: add to second line segment
					GeoLib::Point* intersection (NULL);
					if (node_index_in_ply>0) // test line segment before closest point
						intersection = calcIntersection(line->getNode(0), line->getNode(1), geo_point, ply->getPoint(node_index_in_ply-1));
					if (intersection == NULL && node_index_in_ply<(nLinePnts-1)) // test line segment after closest point
					{
						intersection = calcIntersection(line->getNode(0), line->getNode(1), geo_point, ply->getPoint(node_index_in_ply+1));
						index_offset = 1; // add to second segment
					}
					if (intersection) // intersection found
					{
						if (new_points->size()==319)
							int a=5;
						intersection_count++;
						std::size_t pos = getPointPosInLine((*new_lines)[l], intersection, node_index_in_ply+index_offset-1, line_segment_map[l], eps);
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

	//this->mapOnMesh(mesh);

	this->_geo_objects.addPointVec(new_points, new_geo_name);
	std::vector<size_t> pnt_id_map = this->_geo_objects.getPointVecObj(new_geo_name)->getIDMap();
	for (std::size_t i=0; i<new_lines->size(); ++i)
		(*new_lines)[i]->update(pnt_id_map);
	this->_geo_objects.addPolylineVec(new_lines, new_geo_name);
}

GeoLib::Point* GeoMapper::calcIntersection(GeoLib::Point const*const p1, GeoLib::Point const*const p2, GeoLib::Point const*const q1, GeoLib::Point const*const q2) const
{
	const double x1 = (*p1)[0], x2 = (*p2)[0], x3 = (*q1)[0], x4 = (*q2)[0], z1 = (*p1)[2], z2 = (*p2)[2];
	const double y1 = (*p1)[1], y2 = (*p2)[1], y3 = (*q1)[1], y4 = (*q2)[1];
 
	const double det = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
	if (fabs(det) < std::numeric_limits<double>::epsilon()) return NULL;
 
	const double pre  = (x1*y2 - y1*x2);
	const double post = (x3*y4 - y3*x4);
	const double x = ( pre * (x3 - x4) - (x1 - x2) * post ) / det;
	const double y = ( pre * (y3 - y4) - (y1 - y2) * post ) / det;
 
	int sign = (z1<z2) ? 1:-1;

	// Check if the x and y coordinates are within both line segments
	if (isPntInBoundingBox(x1,y1,x2,y2,x,y) && isPntInBoundingBox(x3,y3,x4,y4,x,y))
	{
		const double denom (fabs(x2-x1));
		if (denom==0) return NULL;
		return new GeoLib::Point(x, y, fabs(fabs(x-x1)*fabs((*p2)[2]-(*p1)[2])/denom) + (*p1)[2]);
	}
	return NULL;
}

std::size_t GeoMapper::getPointPosInLine(GeoLib::Polyline const*const line, GeoLib::Point const*const point, unsigned line_segment, const std::vector<unsigned> &line_segment_map, double eps) const
{
	const std::size_t nPoints (line->getNumberOfPoints());
	const GeoLib::Point* start (line->getPoint(line_segment));
	const double max_dist = MathLib::sqrDist(point, start);
	bool line_segment_found (false);
	for (std::size_t i=0; i<nPoints; ++i)
	{
		if (line_segment_found && line_segment_map[i]>line_segment)
			return i;
		if (line_segment_map[i]==line_segment)
		{
			line_segment_found = true;
			const double v[3] = {(*line->getPoint(i))[0] - (*point)[0], (*line->getPoint(i))[1] - (*point)[1], 0};
			if (MathLib::scpr<double,3>(v,v) < eps) return 0; // don't insert point
			if (MathLib::sqrDist(start, line->getPoint(i))>max_dist)
				return i;
		}
	}
	return 0; // this should not happen!
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

void GeoMapper::mapData(MeshLib::Mesh const*const mesh)
{
	const std::vector<GeoLib::Point*> *points (this->_geo_objects.getPointVec(this->_geo_name));
	GeoLib::Station* stn_test = dynamic_cast<GeoLib::Station*>((*points)[0]);
	bool is_borehole(false);
	if (stn_test != NULL && static_cast<GeoLib::StationBorehole*>((*points)[0])->type() == GeoLib::Station::StationType::BOREHOLE)
		is_borehole = true;
	size_t nPoints (points->size());

	if (!is_borehole)
	{
		for (unsigned j=0; j<nPoints; ++j)
		{
			GeoLib::Point* pnt ((*points)[j]);
			(*pnt)[2] = (_grid) ? this->getMeshElevation((*pnt)[0],(*pnt)[1], mesh)
				                : this->getDemElevation((*pnt)[0],(*pnt)[1]);
		}
	}
	else
	{
		for (unsigned j=0; j<nPoints; ++j)
		{
			GeoLib::Point* pnt ((*points)[j]);
			double offset = (_grid) ? (this->getMeshElevation((*pnt)[0],(*pnt)[1], mesh) - (*pnt)[2])
				                    : (this->getDemElevation((*pnt)[0],(*pnt)[1]) - (*pnt)[2]);

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

float GeoMapper::getDemElevation(double x, double y) const
{
	const double origin_x(_raster->getOrigin()[0]);
	const double origin_y(_raster->getOrigin()[1]);
	const double cellsize(_raster->getRasterPixelDistance());
	const std::size_t width(_raster->getNCols());
	const std::size_t height(_raster->getNRows());

	if ((x<origin_x) || (x>origin_x+(width*cellsize)) || (y<origin_y) || (y>origin_y+(height*cellsize)))
		return 0;

	const unsigned x_index = static_cast<unsigned>((x-origin_x)/cellsize);
	const unsigned y_index = static_cast<unsigned>((y-origin_y)/cellsize);

	return static_cast<float>(*(_raster->begin() + (y_index*width+x_index)));
}

double GeoMapper::getMeshElevation(double x, double y, MeshLib::Mesh const*const mesh) const
{
	double coords[3] = {x,y,0};
	const GeoLib::PointWithID* pnt (_grid->getNearestPoint(coords));
	return (*(mesh->getNode(pnt->getID())))[2];
}

GeoLib::Grid<GeoLib::PointWithID>* GeoMapper::getFlatGrid(MeshLib::Mesh const*const mesh, std::vector<GeoLib::PointWithID*> sfc_pnts) const
{
	if (mesh->getDimension()<3) //much faster
	{
		size_t nNodes (mesh->getNNodes());
		sfc_pnts.resize(nNodes);
		const std::vector<MeshLib::Node*> nodes (mesh->getNodes());
		for (unsigned i(0); i<nNodes; ++i)
			sfc_pnts[i] = new GeoLib::PointWithID(nodes[i]->getCoords(), nodes[i]->getID());
	}
	else
	{
		double dir[3] = {0,0,1};
		sfc_pnts = MeshLib::MeshSurfaceExtraction::getSurfaceNodes(*mesh, dir);
	}
	size_t nPoints (sfc_pnts.size());
	for (unsigned i=0; i<nPoints; ++i)
	{
		GeoLib::PointWithID* pnt (sfc_pnts[i]);
		(*pnt)[2] = 0;
	}

	return new GeoLib::Grid<GeoLib::PointWithID>(sfc_pnts.begin(), sfc_pnts.end());
}

