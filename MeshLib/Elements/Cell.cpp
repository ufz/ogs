/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Implementation of the Cell class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Cell.h"

namespace MeshLib {
/*
Cell::Cell(Node** nodes, MshElemType::type type, unsigned value)
	: Element(nodes, type, value)
{
}
*/
Cell::Cell(unsigned value)
	: Element(value)
{
}

Cell::~Cell()
{}

bool Cell::isOnSurface() const
{
	unsigned n (this->getNNeighbors());
	for (unsigned i(0); i<n; i++)
		if (!this->_neighbors[i])
			return true;
	return false;
}

}

