/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of the MeshSurfaceExtraction class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshSurfaceExtraction.h"

#include <cassert>

#include "logog/include/logog.hpp"

#include "PointWithID.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Face.h"
#include "Elements/Cell.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"

namespace MeshLib {

void MeshSurfaceExtraction::getSurfaceAreaForNodes(const MeshLib::Mesh &mesh, std::vector<double> &node_area_vec)
{
	if (mesh.getDimension() == 2)
	{
		double total_area (0);

		// for each node, a vector containing all the element idget every element
		std::vector<MeshLib::Node*> nodes = mesh.getNodes();
		const size_t nNodes ( mesh.getNNodes() );
		for (size_t n=0; n<nNodes; ++n)
		{
			double node_area (0);

			std::vector<MeshLib::Element*> conn_elems = nodes[n]->getElements();
			const size_t nConnElems (conn_elems.size());

			for (size_t i=0; i<nConnElems; ++i)
			{
				const MeshLib::Element* elem (conn_elems[i]);
				const unsigned nElemParts = (elem->getGeomType() == MeshElemType::TRIANGLE) ? 3 : 4;
				const double area = conn_elems[i]->getContent() / nElemParts;
				node_area += area;
				total_area += area;
			}

			node_area_vec.push_back(node_area);
		}

		INFO ("Total surface Area: %f", total_area);
	}
	else
		ERR ("Error in MeshSurfaceExtraction::getSurfaceAreaForNodes() - Given mesh is no surface mesh (dimension != 2).");
}

MeshLib::Mesh* MeshSurfaceExtraction::getMeshSurface(const MeshLib::Mesh &mesh, const MathLib::Vector3 &dir)
{
	INFO ("Extracting mesh surface...");
	const std::vector<MeshLib::Element*> all_elements (mesh.getElements());
	const std::vector<MeshLib::Node*> all_nodes (mesh.getNodes());

	std::vector<MeshLib::Element*> sfc_elements;
	get2DSurfaceElements(all_elements, sfc_elements, dir, mesh.getDimension());

	if (sfc_elements.empty())
		return nullptr;
	
	std::vector<MeshLib::Node*> sfc_nodes;
	std::vector<unsigned> node_id_map(mesh.getNNodes());
	get2DSurfaceNodes(all_nodes, sfc_nodes, sfc_elements, node_id_map);

	// create new elements vector with newly created nodes
	std::vector<MeshLib::Element*> new_elements;
	new_elements.reserve(sfc_elements.size());
	for (auto elem = sfc_elements.begin(); elem != sfc_elements.end(); ++elem)
	{
		if ((*elem)->getGeomType() == MeshElemType::TRIANGLE) {
			MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
			for (unsigned k(0); k<3; k++)
				tri_nodes[k] = sfc_nodes[node_id_map[(*elem)->getNode(k)->getID()]];
			new_elements.push_back(new MeshLib::Tri(tri_nodes));
		} else {
			assert((*elem)->getGeomType() == MeshElemType::QUAD);
			MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
			for (unsigned k(0); k<4; k++)
				quad_nodes[k] = sfc_nodes[node_id_map[(*elem)->getNode(k)->getID()]];
			new_elements.push_back(new MeshLib::Quad(quad_nodes));
		}
		delete *elem;
	}

	return new Mesh("SurfaceMesh", sfc_nodes, new_elements);
}

void MeshSurfaceExtraction::get2DSurfaceElements(const std::vector<MeshLib::Element*> &all_elements, std::vector<MeshLib::Element*> &sfc_elements, const MathLib::Vector3 &dir, unsigned mesh_dimension)
{
	if (mesh_dimension<2 || mesh_dimension>3)
		ERR("Cannot handle meshes of dimension %i", mesh_dimension);

	bool complete_surface (true);
	if (MathLib::scalarProduct(dir, dir) != 0)
		complete_surface = false;

	for (auto elem = all_elements.begin(); elem != all_elements.end(); ++elem)
	{
		const unsigned element_dimension ((*elem)->getDimension());
		if (element_dimension < mesh_dimension)
			continue;

		if (element_dimension == 2)
		{
			if (!complete_surface)
			{
				MeshLib::Face* face = dynamic_cast<MeshLib::Face*>(*elem);

				if (MathLib::scalarProduct(face->getSurfaceNormal(), dir) <= 0)
					continue;	
			}
			sfc_elements.push_back(*elem);
		}
		else
		{
			const MeshLib::Cell* cell = static_cast<MeshLib::Cell*>(*elem);
			if (!cell->isOnSurface())
				continue;
			const unsigned nFaces (cell->getNFaces());
			for (unsigned j=0; j<nFaces; ++j)
			{
				if (cell->getNeighbor(j) != nullptr)
					continue;

				const MeshLib::Face* face = static_cast<const MeshLib::Face*>(cell->getFace(j));
				if (!complete_surface)
					if (MathLib::scalarProduct(face->getSurfaceNormal(), dir) <= 0)
					{
						delete face;
						continue;
					}
				if (face->getGeomType() == MeshElemType::TRIANGLE)
					sfc_elements.push_back(new MeshLib::Tri(*static_cast<const MeshLib::Tri*>(face)));
				else
					sfc_elements.push_back(new MeshLib::Quad(*static_cast<const MeshLib::Quad*>(face)));
			}
		}
	}
}

void MeshSurfaceExtraction::get2DSurfaceNodes(const std::vector<MeshLib::Node*> &all_nodes, std::vector<MeshLib::Node*> &sfc_nodes, const std::vector<MeshLib::Element*> &sfc_elements, std::vector<unsigned> &node_id_map)
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

std::vector<GeoLib::PointWithID*> MeshSurfaceExtraction::getSurfaceNodes(const MeshLib::Mesh &mesh, const MathLib::Vector3 &dir)
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
