/**
 * \file MshEditor.cpp
 * 2011/06/15 KR Initial implementation
 */

#include "MshEditor.h"
#include "PointWithID.h"
#include "Mesh.h"


MeshLib::Mesh* MshEditor::removeMeshNodes(MeshLib::Mesh* mesh,
                                             const std::vector<size_t> &nodes)
{
	MeshLib::Mesh* new_mesh (new MeshLib::Mesh(*mesh));

	/* TODO6
	// delete nodes and their connected elements and replace them with null pointers
	size_t delNodes = nodes.size();
	for (size_t i = 0; i < delNodes; i++)
	{
		const MeshLib::Node* node = new_mesh->getNode(i);
		std::vector<size_t> conn_elems = node->getConnectedElementIDs();
		for (size_t j = 0; j < conn_elems.size(); j++)
		{
			delete new_mesh->ele_vector[conn_elems[j]];
			new_mesh->ele_vector[conn_elems[j]] = NULL;
		}
		delete new_mesh->nod_vector[nodes[i]];
		new_mesh->nod_vector[nodes[i]] = NULL;
	}

	// create map to adjust node indices in element vector
	size_t nNodes = new_mesh->nod_vector.size();
	std::vector<int> id_map;
	size_t count = 0;
	for (size_t i = 0; i < nNodes; i++)
	{
		if (new_mesh->nod_vector[i])
		{
			new_mesh->nod_vector[i]->SetIndex(count);
			id_map.push_back(count);
			count++;
		}
		else
			id_map.push_back(-1);
	}

	// erase null pointers from node- and element vectors
	for (std::vector<MeshLib::CElem*>::iterator it = new_mesh->ele_vector.begin();
	     it != new_mesh->ele_vector.end(); )
	{
		if (*it)
			++it;
		else
			it = new_mesh->ele_vector.erase(it);
	}

	for (std::vector<MeshLib::CNode*>::iterator it = new_mesh->nod_vector.begin();
	     it != new_mesh->nod_vector.end(); )
	{
		if (*it)
			++it;
		else
			it = new_mesh->nod_vector.erase(it);
	}

	// re-adjust node indices
	size_t nElems = new_mesh->ele_vector.size();
	for (size_t i = 0; i < nElems; i++)
	{
		MeshLib::CElem* elem = new_mesh->ele_vector[i];
		size_t nElemNodes = elem->GetNodesNumber(false);
		for (size_t j = 0; j < nElemNodes; j++)
			elem->SetNodeIndex(j, id_map[elem->GetNodeIndex(j)]);
	}
	*/
	return new_mesh;
}

std::vector<GeoLib::PointWithID*> MshEditor::getSurfaceNodes(const MeshLib::Mesh &mesh)
{
	/* TODO6
	std::cout << "Extracting surface nodes..." << std::endl;
	// Sort points lexicographically
	size_t nNodes (mesh.nod_vector.size());
	std::vector<GeoLib::PointWithID*> nodes;
	std::vector<size_t> perm;
	for (size_t j(0); j<nNodes; j++)
	{
		nodes.push_back(new GeoLib::PointWithID(mesh.nod_vector[j]->getData(), j));
		perm.push_back(j);
	}
	Quicksort<GeoLib::PointWithID*> (nodes, 0, nodes.size(), perm);

	// Extract surface points
	double eps (std::numeric_limits<double>::epsilon());
	*/
	std::vector<GeoLib::PointWithID*> surface_pnts;
	/* TODO6
	for (size_t k(1); k < nNodes; k++)
	{
		const GeoLib::PointWithID& p0 (*(nodes[k - 1]));
		const GeoLib::PointWithID& p1 (*(nodes[k]));
		if (fabs (p0[0] - p1[0]) > eps || fabs (p0[1] - p1[1]) > eps)
			surface_pnts.push_back (nodes[k - 1]);
	}
	// Add last point
	surface_pnts.push_back (nodes[nNodes - 1]);

	*/
	return surface_pnts;
}

MeshLib::Mesh* MshEditor::getMeshSurface(const MeshLib::Mesh &mesh)
{
	/* TODO6
	std::cout << "Extracting mesh surface..." << std::endl;
	GridAdapter surface;
	const std::vector<GeoLib::PointWithID*> sfc_points = MshEditor::getSurfaceNodes(mesh);
	const size_t nSurfacePoints (sfc_points.size());

	std::vector<GridAdapter::Element*> *elements = new std::vector<GridAdapter::Element*>;

	const size_t nElements = mesh.ele_vector.size();
	for (size_t j=0; j<nElements; j++)
	{
		MeshLib::CElem* elem (mesh.ele_vector[j]);
		std::vector<size_t> elem_nodes;
		bool is_surface (true);
		for (size_t i=0; i<4; i++)
		{
			size_t node_index = elem->GetNodeIndex(i);
			bool node_found(false), one_node(true);
			for (size_t k=0; k<nSurfacePoints; k++)
			{
				if (sfc_points[k]->getID() == node_index)
				{
					node_found=true;
					elem_nodes.push_back(k);
					break;
				}
			}
			if (!node_found)
			{
				if (one_node == true)
					one_node = false;
				else
				{
					is_surface = false;
					break;
				}
			}
		}
		if (is_surface)
		{
			GridAdapter::Element* element = new GridAdapter::Element;
			element->material = 0;
			element->type = MshElemType::TRIANGLE;
			element->nodes = elem_nodes;
			elements->push_back(element);
		}
	}

	std::vector<GeoLib::Point*> *nodes = new std::vector<GeoLib::Point*>(nSurfacePoints);
	for (size_t j=0; j<nSurfacePoints; j++)
		//(*nodes)[sfc_points[j]->getID()]=sfc_points[j];
		(*nodes)[j] = sfc_points[j];

	surface.setNodeVector(nodes);
	surface.setElements(elements);

	MeshLib::Mesh* sfc_mesh = new MeshLib::Mesh(*surface.getMesh());
	return sfc_mesh;
	*/
	return new MeshLib::Mesh(mesh);
}

