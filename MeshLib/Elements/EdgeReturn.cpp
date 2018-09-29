/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EdgeReturn.h"

#include <logog/include/logog.hpp>

#include "MeshLib/Node.h"
#include "Element.h"
#include "Line.h"

namespace MeshLib
{

const Element* LinearEdgeReturn::getEdge(const Element* e, unsigned i)
{
    if (i < e->getNumberOfEdges())
    {
        auto** nodes = new Node*[2];
        nodes[0] = const_cast<Node*>(e->getEdgeNode(i,0));
        nodes[1] = const_cast<Node*>(e->getEdgeNode(i,1));
        return new Line(nodes, e->getID());
    }
    ERR("Error in MeshLib::Element::getEdge() - Index does not exist.");
    return nullptr;
}

const Element* QuadraticEdgeReturn::getEdge(const Element* e, unsigned i)
{
    if (i < e->getNumberOfEdges())
    {
        auto** nodes = new Node*[3];
        nodes[0] = const_cast<Node*>(e->getEdgeNode(i,0));
        nodes[1] = const_cast<Node*>(e->getEdgeNode(i,1));
        nodes[2] = const_cast<Node*>(e->getEdgeNode(i,2));
        return new Line3(nodes, e->getID());
    }
    ERR("Error in MeshLib::Element::getEdge() - Index does not exist.");
    return nullptr;
}

} // end MeshLib

