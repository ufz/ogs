/**
 * \file   MeshValidation.cpp
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of the MeshValidation class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshValidation.h"

#include <numeric>

#include "logog/include/logog.hpp"

#include "StringTools.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "MeshEditing/removeMeshNodes.h"
#include "MeshEditing/MeshRevision.h"

namespace MeshLib {

MeshValidation::MeshValidation(MeshLib::Mesh &mesh)
{
	INFO ("Mesh Quality Control:");
	this->removeUnusedMeshNodes(mesh);
	MeshRevision rev(mesh);
	INFO ("Found %d potentially collapsable nodes.", rev.getNCollapsableNodes());

	const std::vector<ElementErrorCode> codes (this->testElementGeometry(mesh));
	std::array<std::string, static_cast<std::size_t>(ElementErrorFlag::MaxValue)> output_str (this->ElementErrorCodeOutput(codes));
	for (std::size_t i = 0; i < output_str.size(); ++i)
		INFO (output_str[i].c_str());
}

std::vector<std::size_t> MeshValidation::findUnusedMeshNodes(const MeshLib::Mesh &mesh)
{
	INFO ("Looking for unused mesh nodes...");
	unsigned count(0);
	const size_t nNodes (mesh.getNNodes());
	const std::vector<MeshLib::Node*> &nodes (mesh.getNodes());
	std::vector<std::size_t> del_node_idx;

	for (unsigned i=0; i<nNodes; ++i)
		if (nodes[i]->getNElements() == 0)
		{
			del_node_idx.push_back(i);
			++count;
		}

	std::string nUnusedNodesStr = (count) ? BaseLib::number2str(count) : "No";
	INFO ("%s unused mesh nodes found.", nUnusedNodesStr.c_str());
	return del_node_idx;
}

std::vector<std::size_t> MeshValidation::removeUnusedMeshNodes(MeshLib::Mesh &mesh)
{
	std::vector<std::size_t> del_node_idx = MeshValidation::findUnusedMeshNodes(mesh);
	MeshLib::removeMeshNodes(mesh, del_node_idx);

	if (!del_node_idx.empty())
		INFO("Removed %d unused mesh nodes.", del_node_idx.size());

	return del_node_idx;
}

 std::vector<ElementErrorCode> MeshValidation::testElementGeometry(const MeshLib::Mesh &mesh)
{
	INFO ("Testing mesh element geometry:");
	const std::size_t nErrorCodes (static_cast<std::size_t>(ElementErrorFlag::MaxValue));
	unsigned error_count[nErrorCodes];
	std::fill_n(error_count, 4, 0);
	const std::size_t nElements (mesh.getNElements());
	const std::vector<MeshLib::Element*> &elements (mesh.getElements());
	std::vector<ElementErrorCode> error_code_vector;
	error_code_vector.reserve(nElements);

	for (std::size_t i=0; i<nElements; ++i)
	{
		const ElementErrorCode e = elements[i]->validate();
		error_code_vector.push_back(e);
		if (e.none())
			continue;

		// increment error statistics
		const std::bitset< static_cast<std::size_t>(ElementErrorFlag::MaxValue) > flags (static_cast< std::bitset<static_cast<std::size_t>(ElementErrorFlag::MaxValue)> >(e));
		for (unsigned j=0; j<nErrorCodes; ++j)
			error_count[j] += flags[j];
	}

	// output
	const unsigned error_sum (static_cast<unsigned>(std::accumulate(error_count, error_count+nErrorCodes, 0.0)));
	if (error_sum != 0)
	{
		ElementErrorFlag flags[nErrorCodes] = { ElementErrorFlag::ZeroVolume, ElementErrorFlag::NonCoplanar, 
											    ElementErrorFlag::NonConvex,  ElementErrorFlag::NodeOrder };
		for (std::size_t i=0; i<nErrorCodes; ++i)
			if (error_count[i])
				INFO ("%d elements found with %s.", error_count[i], ElementErrorCode::toString(flags[i]).c_str());
	}
	else
		INFO ("No errors found.");
	return error_code_vector;
}

std::array<std::string, static_cast<std::size_t>(ElementErrorFlag::MaxValue)>
MeshValidation::ElementErrorCodeOutput(const std::vector<ElementErrorCode> &error_codes)
{
	const std::size_t nErrorFlags (static_cast<std::size_t>(ElementErrorFlag::MaxValue)); 
	ElementErrorFlag flags[nErrorFlags] = { ElementErrorFlag::ZeroVolume, ElementErrorFlag::NonCoplanar, 
		                                    ElementErrorFlag::NonConvex,  ElementErrorFlag::NodeOrder };
	const std::size_t nElements (error_codes.size());
	std::array<std::string, static_cast<std::size_t>(ElementErrorFlag::MaxValue)> output;

	for (std::size_t i=0; i<nErrorFlags; ++i)
	{
		unsigned count(0);
		std::string elementIdStr("");
		
		for (std::size_t j=0; j<nElements; ++j)
		{
			if (error_codes[j][flags[i]])
			{
				elementIdStr += (BaseLib::number2str(j) + ", ");
				count++;
			}
		
		}
		const std::string nErrorsStr = (count) ? BaseLib::number2str(count) : "No";
		output[i] = (nErrorsStr + " elements found with " + ElementErrorCode::toString(flags[i]) + ".\n");

		if (count)
			output[i] += ("ElementIDs: " + elementIdStr + "\n");
	}
	return output;
}

} // end namespace MeshLib
