/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file MshEditor.cpp
 *
 * Created on 2011-06-15 by Karsten Rink
 */

#include "MshEditor.h"
#include "PointWithID.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "Elements/Face.h"
#include "Elements/Cell.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"

#include "MathTools.h"

namespace MeshLib {

void MshEditor::getSurfaceAreaForNodes(const MeshLib::Mesh* mesh, std::vector<double> &node_area_vec)
{
	if (mesh->getDimension() == 2)
	{
		double total_area (0);

		// for each node, a vector containing all the element idget every element
		std::vector<MeshLib::Node*> nodes = mesh->getNodes();
		const size_t nNodes ( mesh->getNNodes() );
		for (size_t n=0; n<nNodes; n++)
		{
			double node_area (0);

			std::vector<MeshLib::Element*> conn_elems = nodes[n]->getElements();
			const size_t nConnElems (conn_elems.size());

			for (size_t i=0; i<nConnElems;i++)
			{
				const MeshLib::Element* elem (conn_elems[i]);
				const unsigned nElemParts = (elem->getType() == MshElemType::TRIANGLE) ? 3 : 4;
				const double area = conn_elems[i]->getContent() / nElemParts;
				node_area += area;
				total_area += area;
			}

			node_area_vec.push_back(node_area);
		}

		std::cout<< "Total surface Area: " << total_area << std::endl;
	}
	else
		std::cout << "Error in MshEditor::getSurfaceAreaForNodes() - Given mesh is no surface mesh (dimension != 2)." << std::endl;
}

MeshLib::Mesh* MshEditor::removeMeshNodes(MeshLib::Mesh* mesh,
                                             const std::vector<size_t> &nodes)
{
	MeshLib::Mesh* new_mesh (new MeshLib::Mesh(*mesh));

	// delete nodes and their connected elements and replace them with null pointers
	const size_t delNodes = nodes.size();
	std::vector<MeshLib::Node*> mesh_nodes = new_mesh->getNodes();
	for (size_t i = 0; i < delNodes; i++)
	{
		const MeshLib::Node* node = new_mesh->getNode(i);
		std::vector<MeshLib::Element*> conn_elems = node->getElements();

		for (size_t j = 0; j < conn_elems.size(); j++)
		{
			delete conn_elems[j];
			conn_elems[j] = NULL;
		}
		delete mesh_nodes[i];
		mesh_nodes[i] = NULL;
	}

	// create map to adjust node indices in element vector
	const size_t nNodes = new_mesh->getNNodes();
	std::vector<int> id_map(nNodes, -1);
	size_t count(0);
	for (size_t i = 0; i < nNodes; i++)
	{
		if (mesh_nodes[i])
		{
			mesh_nodes[i]->setID(count);
			id_map.push_back(count++);
		}
	}

	// erase null pointers from node- and element vectors
	std::vector<MeshLib::Element*> elements = new_mesh->getElements();
	for (std::vector<MeshLib::Element*>::iterator it = elements.begin(); it != elements.end(); )
	{
		if (*it)
			++it;
		else
			it = elements.erase(it);
	}

	for (std::vector<MeshLib::Node*>::iterator it = mesh_nodes.begin(); it != mesh_nodes.end(); )
	{
		if (*it)
			++it;
		else
			it = mesh_nodes.erase(it);
	}

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

MeshLib::Mesh* MshEditor::getMeshSurface(const MeshLib::Mesh &mesh, const double* dir)
{
	std::cout << "Extracting mesh surface..." << std::endl;

	const std::vector<MeshLib::Element*> elements = mesh.getElements();
	std::vector<MeshLib::Element*> new_elements;
	const size_t nElements (mesh.getNElements());

	bool complete_surface = ((dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]) == 0) ? true : false;

	// 2D meshes
	if (mesh.getDimension() == 2 )
	{
		if (complete_surface) return new MeshLib::Mesh(mesh); // simply copy the mesh
		else	// check only surface normal directions of all elements
		{
			for (unsigned i=0; i<nElements; i++)
			{
				MeshLib::Face* face = dynamic_cast<MeshLib::Face*>(elements[i]);
				double normal[3];
				face->getSurfaceNormal(normal);
				if (MathLib::scpr(normal, dir, 3) > 0)
					new_elements.push_back(static_cast<MeshLib::Element*>(face));
			}
		}
	}
	// 3D meshes
	else if (mesh.getDimension() == 3)	//
	{
		for (unsigned i=0; i<nElements; i++)
		{
			if (const MeshLib::Cell* cell = dynamic_cast<MeshLib::Cell*>(elements[i]))
			{
				if (cell->isOnSurface())
				{
					const unsigned nFaces (cell->getNFaces());
					for (unsigned j=0; j<nFaces; j++)
					{
						if (cell->getNeighbor(i) == NULL)
						{
							const MeshLib::Face* face = static_cast<const MeshLib::Face*>(cell->getFace(i));
							if (!complete_surface)
							{
								double normal[3];
								face->getSurfaceNormal(normal);
								if (MathLib::scpr<double,3>(normal, dir) <= 0)
									continue;
							}

							if (face->getType() == MshElemType::TRIANGLE)
								new_elements.push_back(new MeshLib::Tri(*static_cast<const MeshLib::Tri*>(face)));
							else
								new_elements.push_back(new MeshLib::Quad(*static_cast<const MeshLib::Quad*>(face)));
						}
					}
				}
			}
		}

		// now copy nodes
		const size_t nNewElements (new_elements.size());
		std::vector<const MeshLib::Node*> tmp_nodes(mesh.getNNodes(), NULL);
		const size_t nNodes (tmp_nodes.size());
		for (unsigned i=0; i<nNewElements; i++)
		{
			const MeshLib::Element* elem (new_elements[i]);
			for (unsigned j=0; j<elem->getNNodes(); j++)
			{
				const MeshLib::Node* node (elem->getNode(i));
				tmp_nodes[node->getID()] = node;
			}
		}
		std::vector<MeshLib::Node*> new_nodes;
		for (unsigned i=0; i<nNodes; i++)
			if (tmp_nodes[i])
				new_nodes.push_back(new MeshLib::Node(tmp_nodes[i]->getCoords()));
	}


}

} // end namespace MeshLib
