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
#include "Node.h"
#include "Vector3.h"

namespace MeshLib {

const unsigned Cell::dimension = 3u;

/*
Cell::Cell(Node** nodes, MeshElemType type, unsigned value)
	: Element(nodes, type, value)
{
}
*/
Cell::Cell(unsigned value, std::size_t id)
	: Element(value, id), _volume(-1.0) // init with invalid value to detect errors
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

bool Cell::testElementNodeOrder() const
{
	const MathLib::Vector3 c (getCenterOfGravity());
	const unsigned nFaces (this->getNFaces());
	for (unsigned j=0; j<nFaces; ++j)
	{
		MeshLib::Face const*const face (dynamic_cast<const MeshLib::Face*>(this->getFace(j)));
		const MeshLib::Node x (*(face->getNode(1)));
		const MathLib::Vector3 cx (c, x);
		const double s = MathLib::scalarProduct(face->getSurfaceNormal(), cx);
		delete face;
		if (s >= 0)
			return false;
	}
	return true;
}

}

