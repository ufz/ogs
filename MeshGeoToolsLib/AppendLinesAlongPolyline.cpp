
/**
 * @copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "AppendLinesAlongPolyline.h"

#include <logog/include/logog.hpp>

#include "GeoLib/Polyline.h"
#include "GeoLib/PolylineVec.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"

#include "MeshGeoToolsLib/MeshNodesAlongPolyline.h"

namespace MeshGeoToolsLib
{

MeshLib::Mesh* appendLinesAlongPolylines(const MeshLib::Mesh &mesh, const GeoLib::PolylineVec &ply_vec)
{

	// copy existing nodes and elements
	std::vector<MeshLib::Node*> vec_new_nodes = MeshLib::copyNodeVector(mesh.getNodes());
	std::vector<MeshLib::Element*> vec_new_eles = MeshLib::copyElementVector(mesh.getElements(), vec_new_nodes);
	unsigned max_matID = 0;

	auto materialIds = mesh.getProperties().getPropertyVector<int>("MaterialIDs");
	if (materialIds)
	{
		for (auto e : vec_new_eles)
			max_matID = std::max(max_matID, static_cast<unsigned>((*materialIds)[e->getID()]));
	}

	const std::size_t n_ply (ply_vec.size());
	// for each polyline
	for (std::size_t k(0); k < n_ply; k++)
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
		for (std::size_t i=0; i<vec_nodes_on_ply.size()-1; i++) {
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

