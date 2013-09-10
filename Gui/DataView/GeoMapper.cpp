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


void GeoMapper::advancedMapOnMesh(const MeshLib::Mesh* mesh)
{
	const std::vector<GeoLib::Point*> *points (this->_geo_objects.getPointVec(this->_geo_name));
	const std::vector<GeoLib::Polyline*> *plys (this->_geo_objects.getPolylineVec(this->_geo_name));

	const unsigned nGeoPoints ( points->size() );
	std::vector<GeoLib::Point*> *new_points = new std::vector<GeoLib::Point*>(nGeoPoints);
	for (size_t i=0; i<nGeoPoints; ++i)
		(*new_points)[i] = new GeoLib::PointWithID((*(*points)[i])[0],(*(*points)[i])[1],0,i);


	GeoLib::Grid<GeoLib::PointWithID> grid(new_points->begin(), new_points->end());
	const double max_segment_length (this->getMaxSegmentLength(*plys));
	
	const unsigned nMeshNodes ( mesh->getNNodes() );	
	std::vector<int> closest_geo_point(nMeshNodes); 
	std::vector<double> dist(nMeshNodes); 
	for (size_t i=0; i<nMeshNodes; ++i)
	{
		const double* coords = mesh->getNode(i)->getCoords();
		GeoLib::PointWithID* pnt = grid.getNearestPoint(coords);
		dist[i] = MathLib::sqrDist(pnt->getCoords(), coords); //mesh coords z muss noch auf 0 gesetzt werden
		closest_geo_point[i] = (dist[i]<=max_segment_length) ? pnt->getID() : -1;
	}

	std::size_t nLines = plys->size();
	for (std::size_t i=0; i<nMeshNodes; ++i)
	{
		if (closest_geo_point[i] == -1) continue; // is mesh node theoretically close enough?
		const MeshLib::Node* node (mesh->getNode(i));
		if (dist[i] < std::numeric_limits<float>::epsilon()) // is mesh node == geo point?
		{
			(*(*new_points)[closest_geo_point[i]])[2] = (*node)[2];
			continue;
		}

		for (std::size_t l=0; l<nLines; ++l)
		{
			// find relevant polylines
			if (!(*plys)[l]->isPointIDInPolyline(closest_geo_point[i])) continue;
			
			GeoLib::Polyline* ply ((*plys)[l]);
			const std::size_t nLinePoints (ply->getNumberOfPoints());
			std::size_t node_index_in_ply (0);
			for (std::size_t node_index_in_ply=0; node_index_in_ply<nLinePoints; ++node_index_in_ply)
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
					GeoLib::Point* intersection (NULL);
					if (node_index_in_ply>0) // test line segment before closest point
						intersection = calcIntersection(line->getNode(0), line->getNode(1), geo_point, ply->getPoint(node_index_in_ply-1));
					if (intersection == NULL && node_index_in_ply<(nLinePoints-1)) // test line segment after closest point
						intersection = calcIntersection(line->getNode(0), line->getNode(1), geo_point, ply->getPoint(node_index_in_ply+1));
					if (intersection) // intersection found
					{
						intersection_count++;
						new_points->push_back(new GeoLib::PointWithID(intersection->getCoords(), new_points->size()));
						delete intersection;
						//schnittpunkt in linie AN RICHTIGER STELLE einfügen (abstand von neuem punkt zu punkt davor und danach)
					}
				}
			}
		}
	}

	std::string name ("new_points");
	this->_geo_objects.addPointVec(new_points, name);
}


GeoLib::Point* GeoMapper::calcIntersection(GeoLib::Point const*const p1, GeoLib::Point const*const p2, GeoLib::Point const*const q1, GeoLib::Point const*const q2) const
{
	double x1 = (*p1)[0], x2 = (*p2)[0], x3 = (*q1)[0], x4 = (*q2)[0];
	double y1 = (*p1)[1], y2 = (*p2)[1], y3 = (*q1)[1], y4 = (*q2)[1];
 
	double det = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
	if (fabs(det) < std::numeric_limits<double>::epsilon()) return NULL;
 
	// Get the x and y
	double pre  = (x1*y2 - y1*x2);
	double post = (x3*y4 - y3*x4);
	double x = ( pre * (x3 - x4) - (x1 - x2) * post ) / det;
	double y = ( pre * (y3 - y4) - (y1 - y2) * post ) / det;
 
	// Check if the x and y coordinates are within both line segments
	if (isPntInBoundingBox(x1,y1,x2,y2,x,y) && isPntInBoundingBox(x3,y3,x4,y4,x,y))
	{
		double denom (fabs(x2-x1) + (*p1)[2]);
		if (denom==0) return NULL;
		return new GeoLib::Point(x, y, fabs(fabs(x-x1)*fabs((*p2)[2]-(*p1)[2])/denom));
	}
	return NULL;
}

bool GeoMapper::isNodeOnLine(GeoLib::Point const*const p1, GeoLib::Point const*const q1, GeoLib::Point const*const q2) const
{
	double x = (*p1)[0], x1 = (*q1)[0], x2 = (*q2)[0];
	double y = (*p1)[1], y1 = (*q1)[1], y2 = (*q2)[1];

	double det = (x1 - x2) * (y1 - y2) - (y1 - y2) * (x - x2);
	if (fabs(det) < std::numeric_limits<double>::epsilon() && isPntInBoundingBox(x1,y1,x2,y2,x,y))
		return true;
	return false;
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
		double dist(0);
		GeoLib::Polyline* line = lines[i];
		std::size_t nPlyPoints = line->getNumberOfPoints();
		double old_length (0);
		for (size_t j=1; j<nPlyPoints; ++j)
		{
			dist = (line->getLength(j)-old_length);
			old_length = line->getLength(j);
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
