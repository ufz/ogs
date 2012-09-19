/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Face.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Face.h"
#include "Edge.h"
#include "Node.h"

#include "MathTools.h"


namespace MeshLib {
/*
Face::Face(Node** nodes, MshElemType::type type, unsigned value)
	: Element(nodes, type, value)
{
}
*/
Face::Face(unsigned value)
	: Element(value)
{
}

Face::~Face()
{}

void Face::getSurfaceNormal(double normal[3]) const
{
	const double edge1[3] = { (*this->_nodes[0])[0]-(*this->_nodes[1])[0],
				 			  (*this->_nodes[0])[1]-(*this->_nodes[1])[1],
							  (*this->_nodes[0])[2]-(*this->_nodes[1])[2] };
	const double edge2[3] = { (*this->_nodes[1])[0]-(*this->_nodes[2])[0],
							  (*this->_nodes[1])[1]-(*this->_nodes[2])[1],
							  (*this->_nodes[1])[2]-(*this->_nodes[2])[2] };
	MathLib::crossProd(edge1, edge2, normal);
}

}

