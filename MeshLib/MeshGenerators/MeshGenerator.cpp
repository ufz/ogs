/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshGenerator.h"
#include "Node.h"

#include "Elements/Line.h"
#include "Elements/Quad.h"
#include "Elements/Hex.h"

#include <vector>

namespace MeshLib
{
Mesh* MeshGenerator::generateLineMesh(
        const double length,
        const std::size_t subdivision,
        const GeoLib::Point& origin)
{
	const double dx = length / subdivision;

	//nodes
	const std::size_t n_nodes = subdivision + 1;
	std::vector<Node*> nodes(n_nodes);
	for (std::size_t i = 0; i < n_nodes; i++)
		nodes[i] = new Node(origin[0] + dx * i, origin[1], origin[2], i);

	//elements
	const std::size_t n_eles = subdivision;
	std::vector<Element*> elements(n_eles);
	for (std::size_t i = 0; i < n_eles; i++)
	{
		std::array<Node*, 2> element_nodes;
		element_nodes[0] = nodes[i];
		element_nodes[1] = nodes[i + 1];
		elements[i] = new Line(element_nodes);
	}

	return new Mesh("mesh", nodes, elements);
}

Mesh* MeshGenerator::generateRegularQuadMesh(
        const double length,
        const std::size_t subdivision,
        const GeoLib::Point& origin)
{
	return generateRegularQuadMesh(subdivision, subdivision, length/subdivision, origin);
}

Mesh* MeshGenerator::generateRegularQuadMesh(const unsigned n_x_cells,
	                          const unsigned n_y_cells,
                              const double   cell_size,
                              GeoLib::Point const& origin,
							  std::string   const& mesh_name)
{
	//nodes
	const unsigned n_x_nodes (n_x_cells+1);
	const unsigned n_y_nodes (n_y_cells+1);
	std::vector<Node*> nodes;
	nodes.reserve(n_x_nodes * n_y_nodes);
	

	for (std::size_t i = 0, node_id = 0; i < n_y_nodes; i++)
	{
		const double x_offset (origin[0] + cell_size * i);
		for (std::size_t j = 0; j < n_x_nodes; j++)
			nodes.push_back (new Node(origin[1] + cell_size * j, x_offset, origin[2]));
	}

	//elements
	std::vector<Element*> elements;
	elements.reserve(n_x_cells * n_y_cells);

	for (std::size_t j = 0; j < n_y_cells; j++)
	{
		const std::size_t offset_y1 = j * n_x_nodes;
		const std::size_t offset_y2 = (j + 1) * n_x_nodes;
		for (std::size_t k = 0; k < n_x_cells; k++)
		{
			std::array<Node*, 4> element_nodes;
			element_nodes[0] = nodes[offset_y1 + k];
			element_nodes[1] = nodes[offset_y1 + k + 1];
			element_nodes[2] = nodes[offset_y2 + k + 1];
			element_nodes[3] = nodes[offset_y2 + k];
			elements.push_back (new Quad(element_nodes));
		}
	}

	return new Mesh(mesh_name, nodes, elements);
}

Mesh* MeshGenerator::generateRegularHexMesh(
        const double length,
        const std::size_t subdivision,
        const GeoLib::Point& origin)
{
	const double dx = length / subdivision;

	//nodes
	const std::size_t n_nodes_per_dim = subdivision + 1;
	std::vector<Node*> nodes(n_nodes_per_dim * n_nodes_per_dim * n_nodes_per_dim);

	for (std::size_t i = 0, node_id = 0; i < n_nodes_per_dim; i++)
	{
		for (std::size_t j = 0; j < n_nodes_per_dim; j++)
		{
			for (std::size_t k = 0; k < n_nodes_per_dim; k++)
			{
				nodes[node_id] = new Node(origin[0] + dx * k,
										  origin[1] + dx * j,
										  origin[2] + dx * i,
										  node_id);
				node_id++;
			}
		}
	}

	//elements
	const std::size_t n_eles_per_dim = subdivision;
	std::vector<Element*> elements(n_eles_per_dim * n_eles_per_dim * n_eles_per_dim);
	std::size_t elem_id = 0;

	for (std::size_t i = 0; i < n_eles_per_dim; i++)
	{
		const std::size_t offset_z1 = i * n_nodes_per_dim * n_nodes_per_dim; // bottom
		const std::size_t offset_z2 = (i + 1) * n_nodes_per_dim * n_nodes_per_dim; // top
		for (std::size_t j = 0; j < n_eles_per_dim; j++)
		{
			const std::size_t offset_y1 = j * n_nodes_per_dim;
			const std::size_t offset_y2 = (j + 1) * n_nodes_per_dim;
			for (std::size_t k = 0; k < n_eles_per_dim; k++)
			{
				std::array<Node*, 8> element_nodes;
				// bottom
				element_nodes[0] = nodes[offset_z1 + offset_y1 + k];
				element_nodes[1] = nodes[offset_z1 + offset_y1 + k + 1];
				element_nodes[2] = nodes[offset_z1 + offset_y2 + k + 1];
				element_nodes[3] = nodes[offset_z1 + offset_y2 + k];
				// top
				element_nodes[4] = nodes[offset_z2 + offset_y1 + k];
				element_nodes[5] = nodes[offset_z2 + offset_y1 + k + 1];
				element_nodes[6] = nodes[offset_z2 + offset_y2 + k + 1];
				element_nodes[7] = nodes[offset_z2 + offset_y2 + k];
				elements[elem_id++] = new Hex(element_nodes);
			}
		}
	}

	return new Mesh("mesh", nodes, elements);
}

}
