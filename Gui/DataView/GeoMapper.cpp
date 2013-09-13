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
#include "Node.h"
#include "MeshSurfaceExtraction.h"
#include "PointWithID.h"
#include "Raster.h"
#include "readMeshFromFile.h"
#include "StationBorehole.h"


GeoMapper::GeoMapper(GeoLib::GEOObjects &geo_objects, const std::string &geo_name)
	: _geo_objects(geo_objects), _geo_name(geo_name), _mesh(NULL), _grid(NULL), _raster(nullptr)
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

void GeoMapper::mapOnMesh(MeshLib::Mesh* mesh)
{
	this->_mesh = mesh;
	std::vector<GeoLib::PointWithID*> sfc_pnts;
	// init grid
	_grid = this->getFlatGrid(mesh, sfc_pnts);
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
	if (stn_test != NULL && static_cast<GeoLib::StationBorehole*>((*points)[0])->type() == GeoLib::Station::StationType::BOREHOLE)
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
				                : this->getDemElevation((*pnt)[0],(*pnt)[1]);
		}
	}
	else
	{
		for (unsigned j=0; j<nPoints; ++j)
		{
			GeoLib::Point* pnt ((*points)[j]);
			double offset = (_grid) ? (this->getMeshElevation((*pnt)[0],(*pnt)[1], min_val, max_val) - (*pnt)[2])
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

double GeoMapper::getMeshElevation(double x, double y, double min_val, double max_val) const
{
	double coords[3] = {x,y,0};
	const GeoLib::PointWithID* pnt = _grid->getNearestPoint(coords);
	const std::vector<MeshLib::Element*> elements (_mesh->getNode(pnt->getID())->getElements());
	GeoLib::Point* intersection (nullptr);

	for (std::size_t i=0; i<elements.size(); ++i)
	{
		if (intersection==nullptr && elements[i]->getGeomType() == MeshElemType::TRIANGLE)
			intersection=this->triangleLineIntersection(*elements[i]->getNode(0), *elements[i]->getNode(1), *elements[i]->getNode(2), GeoLib::Point(x,y,max_val), GeoLib::Point(x,y,min_val));
	}
	// if the intersection point is not a triangle or something else goes wrong, we simply take the elevation of the nearest point	
	if (intersection)
		return (*intersection)[2];
	return (*(_mesh->getNode(pnt->getID())))[2];
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

GeoLib::Point* GeoMapper::triangleLineIntersection(GeoLib::Point const& a, GeoLib::Point const& b, GeoLib::Point const& c, GeoLib::Point const& p, GeoLib::Point const& q) const
{
	GeoLib::Point pq(q[0]-p[0], q[1]-p[1], q[2]-p[2]);
	GeoLib::Point pa(a[0]-p[0], a[1]-p[1], a[2]-p[2]);
	GeoLib::Point pb(b[0]-p[0], b[1]-p[1], b[2]-p[2]);
	GeoLib::Point pc(c[0]-p[0], c[1]-p[1], c[2]-p[2]);
	
	double u (scalarTriple(pq, pc, pb));
	if (u<0) return nullptr;
	double v (scalarTriple(pq, pa, pc));
	if (v<0) return nullptr;
	double w (scalarTriple(pq, pb, pa));
	if (w<0) return nullptr;
	
	double denom (1.0/(u+v+w));
	u*=denom;
	v*=denom;
	w*=denom;
	return new GeoLib::Point(u*a[0]+v*b[0]+w*c[0],u*a[1]+v*b[1]+w*c[1],u*a[2]+v*b[2]+w*c[2]);
}

double GeoMapper::scalarTriple(GeoLib::Point const& u, GeoLib::Point const& v, GeoLib::Point const& w) const
{
	double cross[3];
	MathLib::crossProd(u.getCoords(), v.getCoords(), cross);
	double result(0);
	for (unsigned i=0; i<3; ++i)
		result+=(cross[i]*w[i]);
	return result;
}


