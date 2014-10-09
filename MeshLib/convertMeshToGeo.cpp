/**
 * \file
 * \author Karsten Rink
 * \date   2013-07-05
 * \brief  Implementation of  of mesh to geometry conversion.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "convertMeshToGeo.h"

#include "logog/include/logog.hpp"

#include "GeoLib/GEOObjects.h"
namespace GeoLib
{
    class Surface;
}

#include "Mesh.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "MeshInformation.h"

namespace MeshLib {

bool convertMeshToGeo(const MeshLib::Mesh &mesh, GeoLib::GEOObjects &geo_objects, double eps)
{
	if (mesh.getDimension() != 2)
	{
		ERR ("Mesh to geometry conversion is only working for 2D meshes.");
		return false;
	}
	
	// nodes to points conversion
	const std::size_t nNodes (mesh.getNNodes());
	std::vector<GeoLib::Point*> *points = new std::vector<GeoLib::Point*>(nNodes);
	const std::vector<MeshLib::Node*> &nodes = mesh.getNodes();

	for (unsigned i=0; i<nNodes; ++i)
		(*points)[i] = new GeoLib::Point(static_cast<GeoLib::Point>(*nodes[i]));

	std::string mesh_name (mesh.getName());
	geo_objects.addPointVec(points, mesh_name, nullptr, eps);
	const std::vector<std::size_t> id_map (geo_objects.getPointVecObj(mesh_name)->getIDMap());

	// elements to surface triangles conversion
	const std::pair<unsigned, unsigned> bounds (MeshInformation::getValueBounds(mesh));
	const unsigned nMatGroups(bounds.second-bounds.first+1);
	std::vector<GeoLib::Surface*> *sfcs = new std::vector<GeoLib::Surface*>;
	sfcs->reserve(nMatGroups);
	for (unsigned i=0; i<nMatGroups; ++i)
		sfcs->push_back(new GeoLib::Surface(*points));

	const std::vector<MeshLib::Element*> &elements = mesh.getElements();
	const std::size_t nElems (mesh.getNElements());

	for (unsigned i=0; i<nElems; ++i)
	{
		MeshLib::Element* e (elements[i]);
		if (e->getGeomType() == MeshElemType::TRIANGLE)
			(*sfcs)[e->getValue()-bounds.first]->addTriangle(id_map[e->getNodeIndex(0)], id_map[e->getNodeIndex(1)], id_map[e->getNodeIndex(2)]);
		if (e->getGeomType() == MeshElemType::QUAD)
		{
			(*sfcs)[e->getValue()-bounds.first]->addTriangle(id_map[e->getNodeIndex(0)], id_map[e->getNodeIndex(1)], id_map[e->getNodeIndex(2)]);
			(*sfcs)[e->getValue()-bounds.first]->addTriangle(id_map[e->getNodeIndex(0)], id_map[e->getNodeIndex(2)], id_map[e->getNodeIndex(3)]);
		}
		// all other element types are ignored (i.e. lines)
	}

	std::for_each(sfcs->begin(), sfcs->end(), [](GeoLib::Surface* sfc){ if (sfc->getNTriangles()==0) delete sfc; sfc = nullptr;});
	auto sfcs_end = std::remove(sfcs->begin(), sfcs->end(), nullptr);
	sfcs->erase(sfcs_end, sfcs->end());

	geo_objects.addSurfaceVec(sfcs, mesh_name);
	return true;
}


}

