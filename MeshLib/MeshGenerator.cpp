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

#include <vector>

namespace MeshLib
{

Mesh* MeshGenerator::generateLineMesh(const double length, const std::size_t subdivision, const double origin_x, const double origin_y, const double origin_z)
{
    const double unit_length = length / subdivision;
    const std::size_t n_nodes_per_axis = subdivision+1;
    const std::size_t n_eles = subdivision;

    //nodes
    std::vector<Node*> nodes;
    std::size_t node_id(0);
    for (std::size_t i_z=0; i_z<n_nodes_per_axis; i_z++) {
        const double x = unit_length*i_z;
        nodes.push_back(new Node(x+origin_x, origin_y, origin_z, node_id++));
    }

    //elements
    std::vector<Element*> elements;
    for (std::size_t i_z=0; i_z<n_eles; i_z++) {
        Node** e_nodes=new Node*[2];
        e_nodes[0] = nodes[i_z];
        e_nodes[1] = nodes[i_z+1];
        elements.push_back(new Edge(e_nodes));
    }

    return new Mesh("mesh", nodes, elements);
}

Mesh* MeshGenerator::generateRegularQuadMesh(const double length, const std::size_t subdivision, const double origin_x, const double origin_y, const double origin_z)
{
    const double unit_length = length / subdivision;
    const std::size_t n_nodes_per_axis = subdivision+1;

    //nodes
    std::vector<Node*> nodes;
    std::size_t node_id(0);
    const double z = origin_z;
    for (std::size_t j_y=0; j_y<n_nodes_per_axis; j_y++) {
        const double y = unit_length*j_y + origin_y;
        for (std::size_t k_x=0; k_x<n_nodes_per_axis; k_x++) {
            const double x = unit_length*k_x + origin_x;
            nodes.push_back(new Node(x, y, z, node_id++));
        }
    }

    //elements
    std::vector<Element*> elements;
    for (std::size_t j=0; j<subdivision; j++) {
        const std::size_t offset_y1 = j*n_nodes_per_axis;
        const std::size_t offset_y2 = (j+1)*n_nodes_per_axis;
        for (std::size_t k=0; k<subdivision; k++) {
            Node** e_nodes=new Node*[4];
            e_nodes[0] = nodes[offset_y1+k];
            e_nodes[1] = nodes[offset_y1+k+1];
            e_nodes[2] = nodes[offset_y2+k+1];
            e_nodes[3] = nodes[offset_y2+k];
            elements.push_back(new Quad(e_nodes));
        }
    }

    return new Mesh("mesh", nodes, elements);
};

}
