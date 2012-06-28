/**
 * \file Face.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Face.h"
#include "Edge.h"

namespace MeshLib {
/*
Face::Face(Node** nodes, MshElemType::type type, unsigned value)
	: Element(nodes, type, value)
{
}
*/
Face::Face(MshElemType::type type, unsigned value)
	: Element(type, value)
{
}

Face::~Face()
{
	delete[] this->_neighbors;
}



}

