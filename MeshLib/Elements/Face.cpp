/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Face class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Face.h"
#include "Line.h"
#include "Node.h"

#include "MathTools.h"
#include "Vector3.h"


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
	const MathLib::Vector3 u (*_nodes[1], *_nodes[0]);
	const MathLib::Vector3 v (*_nodes[1], *_nodes[2]);
	return MathLib::crossProduct(u,v);
}

bool Face::testElementNodeOrder() const
{
	MathLib::Vector3 up_vec (0,0,1);
	return (MathLib::scalarProduct(this->getSurfaceNormal(), up_vec) < 0) ? true : false;
}

}

