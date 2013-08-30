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
	const double dx = length / subdivision;

	//nodes
	const std::size_t n_nodes = subdivision + 1;
	std::vector<Node*> nodes(n_nodes * n_nodes);

	for (std::size_t i = 0, node_id = 0; i < n_nodes; i++)
		for (std::size_t j = 0; j < n_nodes; j++)
		{
			nodes[node_id] = new Node(origin[0] + dx * j,
									  origin[1] + dx * i,
									  origin[2],
									  node_id);
			node_id++;
		}

	//elements
	std::size_t const n_eles = subdivision;
	std::vector<Element*> elements(n_eles * n_eles);
	std::size_t elem_id = 0;

	for (std::size_t j = 0; j < subdivision; j++)
	{
		const std::size_t offset_y1 = j * n_nodes;
		const std::size_t offset_y2 = (j + 1) * n_nodes;
		for (std::size_t k = 0; k < subdivision; k++)
		{
			std::array<Node*, 4> element_nodes;
			element_nodes[0] = nodes[offset_y1 + k];
			element_nodes[1] = nodes[offset_y1 + k + 1];
			element_nodes[2] = nodes[offset_y2 + k + 1];
			element_nodes[3] = nodes[offset_y2 + k];
			elements[elem_id++] = new Quad(element_nodes);
		}
	}

	return new Mesh("mesh", nodes, elements);
}
}
