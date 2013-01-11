/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Node class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Node.h"
#include "Elements/Element.h"

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

void Node::updateCoordinates(double x, double y, double z)
{
	_x[0] = x;
	_x[1] = y;
	_x[2] = z;

	const size_t nElements (this->_elements.size());
	for (unsigned i=0; i<nElements; i++)
		_elements[i]->computeVolume();
}

}

