/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file Tet10.cpp
 *
 * Created on 2012-05-03 by Karsten Rink
 */

#include "Tet10.h"
#include "Node.h"

namespace MeshLib {

Tet10::Tet10(Node* nodes[10], unsigned value)
	: Tet(value), FemElem()
{
	_nodes = nodes;
	this->_volume = this->computeVolume();
	this->calcCentroid();
}

Tet10::Tet10(const Tet &tet)
	: Tet(tet.getValue()), FemElem()
{
	_nodes = new Node*[10];
	unsigned nNodes (tet.getNNodes());
	for (unsigned i=0; i<nNodes; i++)
		_nodes[i] = const_cast<Node*>(tet.getNode(i));

	if (nNodes < this->getNNodes())
	{
		//TODO: Calculate additional nodes!
	}

	this->_volume = this->computeVolume();
	this->calcCentroid();
}

Tet10::Tet10(const Tet10 &tet)
	: Tet(tet.getValue()), FemElem()
{
	_nodes = new Node*[10];
	for (unsigned i=0; i<10; i++)
		_nodes[i] = tet._nodes[i];
	_centroid = tet.getCentroid();
}

Tet10::~Tet10()
{
}

void Tet10::calcCentroid()
{
	//TODO calculation!
	_centroid = 0;
}

}

