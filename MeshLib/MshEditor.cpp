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

#include "logog.hpp"

namespace MeshLib {

void MshEditor::getSurfaceAreaForNodes(const MeshLib::Mesh* mesh, std::vector<double> &node_area_vec)
{
	if (mesh->getDimension() == 2)
	{
		double total_area (0);

		// for each node, a vector containing all the element idget every element
		std::vector<MeshLib::Node*> nodes = mesh->getNodes();
		const size_t nNodes ( mesh->getNNodes() );
		for (size_t n=0; n<nNodes; ++n)
		{
			double node_area (0);

			std::vector<MeshLib::Element*> conn_elems = nodes[n]->getElements();
			const size_t nConnElems (conn_elems.size());

			for (size_t i=0; i<nConnElems; ++i)
			{
				const MeshLib::Element* elem (conn_elems[i]);
				const unsigned nElemParts = (elem->getType() == MshElemType::TRIANGLE) ? 3 : 4;
				const double area = conn_elems[i]->getContent() / nElemParts;
				node_area += area;
				total_area += area;
			}

			node_area_vec.push_back(node_area);
		}

		INFO ("Total surface Area: %f", total_area);
	}
	else
		ERR ("Error in MshEditor::getSurfaceAreaForNodes() - Given mesh is no surface mesh (dimension != 2).");
}

MeshLib::Mesh* MshEditor::removeMeshNodes(MeshLib::Mesh* mesh,
                                             const std::vector<size_t> &nodes)
{
	MeshLib::Mesh* new_mesh (new MeshLib::Mesh(*mesh));

	// delete nodes and their connected elements and replace them with null pointers
	const size_t delNodes = nodes.size();
	std::vector<MeshLib::Node*> mesh_nodes = new_mesh->getNodes();
	for (size_t i = 0; i < delNodes; ++i)
	{
		const MeshLib::Node* node = new_mesh->getNode(i);
		std::vector<MeshLib::Element*> conn_elems = node->getElements();

		for (size_t j = 0; j < conn_elems.size(); ++j)
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
	for (size_t i = 0; i < nNodes; ++i)
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

MeshLib::Mesh* MshEditor::getMeshSurface(const MeshLib::Mesh &mesh, const double* dir)
{
	INFO ("Extracting mesh surface...");
	const std::vector<MeshLib::Element*> all_elements (mesh.getElements());
	const std::vector<MeshLib::Node*> all_nodes (mesh.getNodes());

	std::vector<MeshLib::Element*> sfc_elements;
	get2DSurfaceElements(all_elements, sfc_elements, dir, mesh.getDimension());

	std::vector<MeshLib::Node*> sfc_nodes;
	std::vector<unsigned> node_id_map(mesh.getNNodes());
	get2DSurfaceNodes(all_nodes, sfc_nodes, sfc_elements, node_id_map);

	// create new elements vector with newly created nodes
	const size_t nNewElements (sfc_elements.size());
	std::vector<MeshLib::Element*> new_elements(sfc_elements.size());
	for (unsigned i=0; i<nNewElements; ++i)
	{
		MeshLib::Element* elem (sfc_elements[i]);
		if (elem->getType() == MshElemType::TRIANGLE) {
			MeshLib::Node** tri_nodes (new MeshLib::Node*[3]);
			for (unsigned k(0); k<3; k++)
				tri_nodes[k] = sfc_nodes[node_id_map[elem->getNode(k)->getID()]];
			new_elements[i] = new MeshLib::Tri(tri_nodes);
		} else {
			new_elements[i] = new MeshLib::Quad(sfc_nodes[node_id_map[elem->getNode(0)->getID()]],
				                                sfc_nodes[node_id_map[elem->getNode(1)->getID()]],
												sfc_nodes[node_id_map[elem->getNode(2)->getID()]],
												sfc_nodes[node_id_map[elem->getNode(3)->getID()]]);
		}
		delete sfc_elements[i];
	}

	return new Mesh("SurfaceMesh", sfc_nodes, new_elements);
}

void MshEditor::get2DSurfaceElements(const std::vector<MeshLib::Element*> &all_elements, std::vector<MeshLib::Element*> &sfc_elements, const double* dir, unsigned mesh_dimension)
{
	bool complete_surface (true);
	if (dir)
		complete_surface = ((dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]) == 0) ? true : false;

	const size_t nElements (all_elements.size());

	if (mesh_dimension == 2 )
	{
		for (unsigned i=0; i<nElements; ++i)
		{
			if (complete_surface)
				sfc_elements.push_back(all_elements[i]);
			else
			{
				MeshLib::Face* face = dynamic_cast<MeshLib::Face*>(all_elements[i]);
				double normal[3];
				face->getSurfaceNormal(normal);
				if (MathLib::scpr(normal, dir, 3) > 0)
					sfc_elements.push_back(static_cast<MeshLib::Element*>(face));
			}
		}
	}
	else if (mesh_dimension == 3)
	{
		for (unsigned i=0; i<nElements; ++i)
		{
			if (all_elements[i]->getDimension()==3)
			{
				const MeshLib::Cell* cell = static_cast<MeshLib::Cell*>(all_elements[i]);
				if (cell->isOnSurface())
				{
					const unsigned nFaces (cell->getNFaces());
					for (unsigned j=0; j<nFaces; ++j)
					{
						if (cell->getNeighbor(j) == NULL)
						{
							const MeshLib::Face* face = static_cast<const MeshLib::Face*>(cell->getFace(j));
							if (!complete_surface)
							{
								double normal[3];
								face->getSurfaceNormal(normal);
								if (MathLib::scpr<double,3>(normal, dir) <= 0)
									continue;
							}

							if (face->getType() == MshElemType::TRIANGLE)
								sfc_elements.push_back(new MeshLib::Tri(*static_cast<const MeshLib::Tri*>(face)));
							else
								sfc_elements.push_back(new MeshLib::Quad(*static_cast<const MeshLib::Quad*>(face)));
						}
					}
				}
			}
		}
	}
}

void MshEditor::get2DSurfaceNodes(const std::vector<MeshLib::Node*> &all_nodes, std::vector<MeshLib::Node*> &sfc_nodes, const std::vector<MeshLib::Element*> &sfc_elements, std::vector<unsigned> &node_id_map)
{
	const size_t nNewElements (sfc_elements.size());
	std::vector<const MeshLib::Node*> tmp_nodes(all_nodes.size(), NULL);
	const size_t nNodes (tmp_nodes.size());
	for (unsigned i=0; i<nNewElements; ++i)
	{
		const MeshLib::Element* elem (sfc_elements[i]);
		for (unsigned j=0; j<elem->getNNodes(); ++j)
		{
			const MeshLib::Node* node (elem->getNode(j));
			tmp_nodes[node->getID()] = node;
		}
	}
	for (unsigned i=0; i<nNodes; ++i)
	{
		if (tmp_nodes[i])
		{
			node_id_map[i] = sfc_nodes.size();
			sfc_nodes.push_back(new MeshLib::Node(tmp_nodes[i]->getCoords(), tmp_nodes[i]->getID()));
		}
	}
}

std::vector<GeoLib::PointWithID*> MshEditor::getSurfaceNodes(const MeshLib::Mesh &mesh, const double *dir)
{
	INFO ("Extracting surface nodes...");
	const std::vector<MeshLib::Element*> all_elements (mesh.getElements());
	const std::vector<MeshLib::Node*> all_nodes (mesh.getNodes());

	std::vector<MeshLib::Element*> sfc_elements;
	get2DSurfaceElements(all_elements, sfc_elements, dir, mesh.getDimension());

	std::vector<MeshLib::Node*> sfc_nodes;
	std::vector<unsigned> node_id_map(mesh.getNNodes());
	get2DSurfaceNodes(all_nodes, sfc_nodes, sfc_elements, node_id_map);

	const unsigned nElements (sfc_elements.size());
	for (unsigned i=0; i<nElements; ++i)
		delete sfc_elements[i];

	const size_t nNodes (sfc_nodes.size());
	std::vector<GeoLib::PointWithID*> surface_pnts(nNodes);
	for (unsigned i=0; i<nNodes; ++i)
	{
		surface_pnts[i] = new GeoLib::PointWithID(sfc_nodes[i]->getCoords(), sfc_nodes[i]->getID());
		delete sfc_nodes[i];
	}
	return surface_pnts;
}


} // end namespace MeshLib
