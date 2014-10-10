/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshGenerator.h"

#include <vector>

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Tri.h"


namespace MeshLib
{
Mesh* MeshGenerator::generateLineMesh(
        const double length,
        const std::size_t subdivision,
        const GeoLib::Point& origin)
{
	return MeshGenerator::generateLineMesh(subdivision, length/subdivision, origin);
}

Mesh* MeshGenerator::generateLineMesh(const unsigned n_cells,
                       const double   cell_size,
                       GeoLib::Point const& origin,
					   std::string   const& mesh_name)
{
	//nodes
	const std::size_t n_nodes (n_cells + 1);
	std::vector<Node*> nodes;
	nodes.reserve(n_nodes);

	for (std::size_t i = 0; i < n_nodes; i++)
		nodes.push_back (new Node(origin[0] + cell_size * i, origin[1], origin[2], i));

	//elements
	std::vector<Element*> elements;
	elements.reserve(n_cells);

	for (std::size_t i = 0; i < n_cells; i++)
	{
		std::array<Node*, 2> element_nodes;
		element_nodes[0] = nodes[i];
		element_nodes[1] = nodes[i + 1];
		elements.push_back (new Line(element_nodes));
	}

	return new Mesh(mesh_name, nodes, elements);
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
                              const double cell_size,
                              GeoLib::Point const& origin,
                              std::string const& mesh_name)
{
	//nodes
	const unsigned n_x_nodes (n_x_cells+1);
	const unsigned n_y_nodes (n_y_cells+1);
	std::vector<Node*> nodes;
	nodes.reserve(n_x_nodes * n_y_nodes);	

	for (std::size_t i = 0; i < n_y_nodes; i++)
	{
		const double y_offset (origin[1] + cell_size * i);
		for (std::size_t j = 0; j < n_x_nodes; j++)
			nodes.push_back (new Node(origin[0] + cell_size * j, y_offset, origin[2]));
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
	return MeshGenerator::generateRegularHexMesh(subdivision, subdivision, subdivision, length/subdivision, origin);
}

Mesh* MeshGenerator::generateRegularHexMesh(const unsigned n_x_cells,
	                         const unsigned n_y_cells,
	                         const unsigned n_z_cells,
	                         const double   cell_size,
	                         GeoLib::Point const& origin,
	                         std::string   const& mesh_name)
{
	//nodes
	const unsigned n_x_nodes (n_x_cells+1);
	const unsigned n_y_nodes (n_y_cells+1);
	const unsigned n_z_nodes (n_z_cells+1);
	std::vector<Node*> nodes;
	nodes.reserve(n_x_nodes * n_y_nodes * n_z_nodes);

	for (std::size_t i = 0; i < n_z_nodes; i++)
	{
		const double z_offset (origin[2] + cell_size * i);
		for (std::size_t j = 0; j < n_y_nodes; j++)
		{
			const double y_offset (origin[1] + cell_size * j);
			for (std::size_t k = 0; k < n_x_nodes; k++)
				nodes.push_back (new Node(origin[0] + cell_size * k, y_offset, z_offset));
		}
	}

	//elements
	std::vector<Element*> elements;
	elements.reserve(n_x_cells * n_y_cells * n_z_cells);

	for (std::size_t i = 0; i < n_z_cells; i++)
	{
		const std::size_t offset_z1 = i * n_x_nodes * n_y_nodes; // bottom
		const std::size_t offset_z2 = (i + 1) * n_x_nodes * n_y_nodes; // top
		for (std::size_t j = 0; j < n_y_cells; j++)
		{
			const std::size_t offset_y1 = j * n_x_nodes;
			const std::size_t offset_y2 = (j + 1) * n_x_nodes;
			for (std::size_t k = 0; k < n_x_cells; k++)
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
				elements.push_back (new Hex(element_nodes));
			}
		}
	}

	return new Mesh(mesh_name, nodes, elements);
}

Mesh* MeshGenerator::generateRegularTriMesh(
        const double length,
        const std::size_t subdivision,
        const GeoLib::Point& origin)
{
	return generateRegularTriMesh(subdivision, subdivision, length/subdivision, origin);
}

Mesh* MeshGenerator::generateRegularTriMesh(const unsigned n_x_cells,
                              const unsigned n_y_cells,
                              const double cell_size,
                              GeoLib::Point const& origin,
                              std::string const& mesh_name)
{
	//nodes
	const unsigned n_x_nodes (n_x_cells+1);
	const unsigned n_y_nodes (n_y_cells+1);
	std::vector<Node*> nodes;
	nodes.reserve(n_x_nodes * n_y_nodes);

	for (std::size_t i = 0; i < n_y_nodes; i++)
	{
		const double y_offset (origin[1] + cell_size * i);
		for (std::size_t j = 0; j < n_x_nodes; j++)
			nodes.push_back (new Node(origin[0] + cell_size * j, y_offset, origin[2]));
	}

	//elements
	std::vector<Element*> elements;
	elements.reserve(n_x_cells * n_y_cells * 2);

	for (std::size_t j = 0; j < n_y_cells; j++)
	{
		const std::size_t offset_y1 = j * n_x_nodes;
		const std::size_t offset_y2 = (j + 1) * n_x_nodes;
		for (std::size_t k = 0; k < n_x_cells; k++)
		{
			std::array<Node*, 3> element1_nodes;
			element1_nodes[0] = nodes[offset_y1 + k];
			element1_nodes[1] = nodes[offset_y2 + k + 1];
			element1_nodes[2] = nodes[offset_y2 + k];
			elements.push_back (new Tri(element1_nodes));
			std::array<Node*, 3> element2_nodes;
			element2_nodes[0] = nodes[offset_y1 + k];
			element2_nodes[1] = nodes[offset_y1 + k + 1];
			element2_nodes[2] = nodes[offset_y2 + k + 1];
			elements.push_back (new Tri(element2_nodes));
		}
	}

	return new Mesh(mesh_name, nodes, elements);
}

}
