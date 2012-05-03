/**
 * Face.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Face.h"

namespace MeshLib {

Face::Face(Node** nodes, MshElemType::type type, size_t value)
	: Element(nodes, type, value)
{
}

Face::Face(MshElemType::type type, size_t value)
	: Element(type, value)
{
}

Face::~Face()
{
}


}

