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


#include <numeric>

#include "MeshQualityController.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "MeshEditing/removeMeshNodes.h"


#include "logog/include/logog.hpp"

namespace MeshLib {

MeshQualityController::MeshQualityController(MeshLib::Mesh &mesh)
{
	this->removeUnusedMeshNodes(mesh);
	this->testElementGeometry(mesh);
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

void MeshQualityController::testElementGeometry(const MeshLib::Mesh &mesh)
{
	const std::size_t nErrorCodes (static_cast<std::size_t>(ElementErrorFlag::MaxValue));
	unsigned error_count[nErrorCodes] = {{0}};
	const std::size_t nElements (mesh.getNElements());
	const std::vector<MeshLib::Element*> &elements (mesh.getElements());

	unsigned count(0);
	for (unsigned i=0; i<nElements; ++i)
	{
		const ElementErrorCode e = elements[i]->isValid();
		if (e.none())
			continue;

		const std::bitset<nErrorCodes> flags (e.bitset());
		for (unsigned i=0; i<nErrorCodes; ++i)
			error_count[i] += flags[i];
	}

	const unsigned error_sum (static_cast<unsigned>(std::accumulate(error_count, error_count+nErrorCodes, 0.0)));
	if (error_sum != 0)
	{
		ElementErrorFlag flags[nErrorCodes] = {ElementErrorFlag::ZeroVolume, ElementErrorFlag::NonCoplanar, 
											   ElementErrorFlag::NonConvex,  ElementErrorFlag::NodeOrder };
		for (std::size_t i=0; i<nErrorCodes; ++i)
			if (error_count[i])
				INFO ("%d elements found with %s.", error_count[i], ElementErrorCode::toString(flags[i]));
	}
	else
		INFO ("No errors found.");
}

} // end namespace MeshLib
