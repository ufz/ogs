/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Node class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshLib/Node.h"
#include "Elements/Element.h"

namespace MeshLib {

Node::Node(const double coords[3], std::size_t id)
    : MathLib::Point3dWithID(
        std::array<double,3>{{coords[0], coords[1], coords[2]}}, id)
{
}

Node::Node(std::array<double, 3> const& coords, std::size_t id)
    : MathLib::Point3dWithID(coords, id)
{
}

Node::Node(double x, double y, double z, std::size_t id)
    : MathLib::Point3dWithID(std::array<double,3>({{x, y, z}}), id)
{
}

Node::Node(const Node &node)
    : MathLib::Point3dWithID(node.x_, node.getID())
{
}

void Node::updateCoordinates(double x, double y, double z)
{
    x_[0] = x;
    x_[1] = y;
    x_[2] = z;

    const std::size_t nElements (this->elements_.size());
    for (std::size_t i = 0; i < nElements; i++)
    {
        elements_[i]->computeVolume();
    }
}

bool isBaseNode(Node const& node)
{
    // Check if node is connected.
    if (node.getNumberOfElements() == 0)
    {
        return true;
    }

    // In a mesh a node always belongs to at least one element.
    auto const e = node.getElement(0);

    auto const n_base_nodes = e->getNumberOfBaseNodes();
    auto const local_index = e->getNodeIDinElement(&node);
    return local_index < n_base_nodes;
}
}  // namespace MeshLib
