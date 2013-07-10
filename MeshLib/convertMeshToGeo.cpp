/**
 * \file
 * \author Karsten Rink
 * \date   2013-07-05
 * \brief  Implementation of  of mesh to geometry conversion.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "convertMeshToGeo.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "GEOObjects.h"
#include "Surface.h"

#include "Mesh.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Node.h"

namespace MeshLib {

bool convertMeshToGeo(const MeshLib::Mesh &mesh, GeoLib::GEOObjects* geo_objects)
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

	// elements to surface triangles conversion
	const std::vector<MeshLib::Element*> &elements = mesh.getElements();
	GeoLib::Surface* sfc = new GeoLib::Surface(*points);
	const std::size_t nElems (mesh.getNElements());

	for (unsigned i=0; i<nElems; ++i)
	{
		MeshLib::Element* e (elements[i]);
		if (e->getGeomType() == MshElemType::TRIANGLE)
			sfc->addTriangle(e->getNodeIndex(0), e->getNodeIndex(1), e->getNodeIndex(2));
		if (e->getGeomType() == MshElemType::QUAD)
		{
			sfc->addTriangle(e->getNodeIndex(0), e->getNodeIndex(1), e->getNodeIndex(2));
			sfc->addTriangle(e->getNodeIndex(0), e->getNodeIndex(2), e->getNodeIndex(3));
		}
		// all other element types are ignored (i.e. lines)
	}

	std::vector<GeoLib::Surface*> *sfcs = new std::vector<GeoLib::Surface*>(1);
	(*sfcs)[0] = sfc;

	std::string mesh_name (mesh.getName());
	geo_objects->addPointVec(points, mesh_name);
	geo_objects->addSurfaceVec(sfcs, mesh_name);
	return true;
}


}

