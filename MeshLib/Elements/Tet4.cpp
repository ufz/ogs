/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file Tet4.cpp
 *
 * Created on 2012-05-03 by Karsten Rink
 */

#include "Tet4.h"
#include "Node.h"

namespace MeshLib {

Tet4::Tet4(Node* nodes[4], unsigned value)
	: Tet(nodes, value), FemElem()
{
	this->calcCentroid();
}

Tet4::Tet4(const Tet &tet)
	: Tet(tet), FemElem()
{
	this->calcCentroid();
}

Tet4::Tet4(const Tet4 &tet)
	: Tet(tet.getValue()), FemElem()
{
	for (unsigned i=0; i<4; i++)
		_nodes[i] = tet._nodes[i];
	_centroid = tet.getCentroid();
}

Tet4::~Tet4()
{
}

void Tet4::calcCentroid()
{
	_centroid = 0;
	//TODO calculation!
}

}

