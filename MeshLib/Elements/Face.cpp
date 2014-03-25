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
Face::Face(unsigned value, std::size_t id)
	: Element(value, id), _area(-1.0) // init with invalid value to detect errors
{
}

Face::~Face()
{}

MathLib::Vector3 Face::getSurfaceNormal() const
{
	const MathLib::Vector3 u ((*this->_nodes[0])[0]-(*this->_nodes[1])[0],
				 		      (*this->_nodes[0])[1]-(*this->_nodes[1])[1],
						      (*this->_nodes[0])[2]-(*this->_nodes[1])[2]);
	const MathLib::Vector3 v ((*this->_nodes[1])[0]-(*this->_nodes[2])[0],
							  (*this->_nodes[1])[1]-(*this->_nodes[2])[1],
							  (*this->_nodes[1])[2]-(*this->_nodes[2])[2]);
	return MathLib::crossProduct(u,v);
}

}

