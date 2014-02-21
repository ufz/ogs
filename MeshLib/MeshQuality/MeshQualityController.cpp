/**
 * \file   MeshQualityController.cpp
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of the MeshQualityController class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include "MeshQualityController.h"
#include "Mesh.h"
#include "Node.h"
#include "MeshEditing/removeMeshNodes.h"


#include "logog/include/logog.hpp"

namespace MeshLib {

MeshQualityController::MeshQualityController(MeshLib::Mesh &mesh)
{
	this->removeUnusedMeshNodes(mesh);
}

void MeshQualityController::removeUnusedMeshNodes(MeshLib::Mesh &mesh)
{
	std::vector<MeshLib::Node*> nodes (mesh.getNodes());
	std::vector<std::size_t> del_node_idx;
	std::size_t nNodes (mesh.getNNodes());
	for (std::size_t i=0; i<nNodes; ++i)
	{
		if (nodes[i]->getNElements() == 0)
			del_node_idx.push_back(i);
	}
	MeshLib::removeMeshNodes(mesh, del_node_idx);

	if (!del_node_idx.empty())
		INFO("Removed %d unused mesh nodes.", del_node_idx.size());
}

void MeshQualityController::testElementGeometry(MeshLib::Mesh &mesh)
{
	unsigned count(0);
	const std::size_t nElements (mesh.getNElements());
	for (std::size_t i=0; i<nElements; ++i)
	{

	}
}

} // end namespace MeshLib
