
/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "AppendLinesAlongPolyline.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "Elements/Line.h"
#include "Elements/Element.h"
#include "MeshEnums.h"
#include "MeshEditing/DuplicateMeshComponents.h"

#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"

namespace MeshGeoToolsLib
{

MeshLib::Mesh* appendLinesAlongPolylines(const MeshLib::Mesh &mesh, const GeoLib::PolylineVec &ply_vec)
{

	// copy existing nodes and elements
	std::vector<MeshLib::Node*> vec_new_nodes = MeshLib::copyNodeVector(mesh.getNodes());
	std::vector<MeshLib::Element*> vec_new_eles = MeshLib::copyElementVector(mesh.getElements(), vec_new_nodes);
	unsigned max_matID = 0;
	for (auto e : vec_new_eles)
		max_matID = std::max(max_matID, e->getValue());

	const size_t n_ply (ply_vec.size());
	// for each polyline
	for (size_t k(0); k < n_ply; k++)
	{
		const GeoLib::Polyline* ply = (*ply_vec.getVector())[k];

		// search nodes on the polyline
		MeshGeoToolsLib::MeshNodesAlongPolyline mshNodesAlongPoly(mesh, *ply, mesh.getMinEdgeLength()*0.5);
		auto &vec_nodes_on_ply = mshNodesAlongPoly.getNodeIDs();
		if (vec_nodes_on_ply.empty()) {
			std::string ply_name;
			ply_vec.getNameOfElementByID(k, ply_name);
			INFO("No nodes found on polyline %s", ply_name.c_str());
			continue;
		}

		// add line elements
		for (size_t i=0; i<vec_nodes_on_ply.size()-1; i++) {
			std::array<MeshLib::Node*, 2> element_nodes;
			element_nodes[0] = vec_new_nodes[vec_nodes_on_ply[i]];
			element_nodes[1] = vec_new_nodes[vec_nodes_on_ply[i+1]];
			vec_new_eles.push_back(new MeshLib::Line(element_nodes, max_matID+k+1));
		}
	}

	// generate a mesh
	const std::string new_mesh_name = mesh.getName() + "_with_lines";
	return new MeshLib::Mesh(new_mesh_name, vec_new_nodes, vec_new_eles);
}

} // MeshGeoToolsLib

