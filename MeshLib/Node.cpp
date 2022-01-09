/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Node class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshLib/Node.h"

#include "Elements/Element.h"

namespace MeshLib
{
Node::Node(const double coords[3], std::size_t id)
    : MathLib::Point3dWithID(
          std::array<double, 3>{{coords[0], coords[1], coords[2]}}, id)
{
}

Node::Node(std::array<double, 3> const& coords, std::size_t id)
    : MathLib::Point3dWithID(coords, id)
{
}

Node::Node(double x, double y, double z, std::size_t id)
    : MathLib::Point3dWithID(std::array<double, 3>({{x, y, z}}), id)
{
}

Node::Node(const Node& node) : MathLib::Point3dWithID(node, node.getID()) {}

void Node::updateCoordinates(double x, double y, double z)
{
    (*this)[0] = x;
    (*this)[1] = y;
    (*this)[2] = z;
}
}  // namespace MeshLib
