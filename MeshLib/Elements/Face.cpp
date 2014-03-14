/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Face class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Face.h"
#include "Line.h"
#include "Node.h"

#include "MathTools.h"


namespace MeshLib {

const unsigned Face::dimension = 2u;

/*
Face::Face(Node** nodes, MeshElemType type, unsigned value)
	: Element(nodes, type, value)
{
}
*/
Face::Face(unsigned value, unsigned id)
	: Element(value, id), _area(-1.0) // init with invalid value to detect errors
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

