/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file Node.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Node.h"

namespace MeshLib {

Node::Node(const double coords[3], unsigned id)
	: GeoLib::PointWithID(coords, id)
{
}

Node::Node(double x, double y, double z, unsigned id)
	: GeoLib::PointWithID(x, y, z, id)
{
}

Node::Node(const Node &node)
	: GeoLib::PointWithID(node.getCoords(), node.getID())
{
}

Node::~Node()
{
}

}

