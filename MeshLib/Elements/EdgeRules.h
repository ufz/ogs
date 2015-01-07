/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EDGE_RULES_H_
#define EDGE_RULES_H_

#include "logog/include/logog.hpp"

#include "MeshLib/Node.h"
#include "Element.h"
#include "Line.h"

namespace MeshLib
{

namespace detail
{

class LinearEdgeReturn
{
public:
	static const Element* getEdge(const Element* e, unsigned i)
	{
		if (i < e->getNEdges())
		{
			Node** nodes = new Node*[2];
			nodes[0] = const_cast<Node*>(e->getEdgeNode(i,0));
			nodes[1] = const_cast<Node*>(e->getEdgeNode(i,1));
			return new Line(nodes);
		}
		ERR("Error in MeshLib::Element::getEdge() - Index does not exist.");
		return nullptr;
	}
};

class QuadraticEdgeReturn
{
public:
	static const Element* getEdge(const Element* e, unsigned i)
	{
		if (i < e->getNEdges())
		{
			Node** nodes = new Node*[3];
			nodes[0] = const_cast<Node*>(e->getEdgeNode(i,0));
			nodes[1] = const_cast<Node*>(e->getEdgeNode(i,1));
			nodes[2] = const_cast<Node*>(e->getEdgeNode(i,2));
			return new Line3(nodes);
		}
		ERR("Error in MeshLib::Element::getEdge() - Index does not exist.");
		return nullptr;
	}
};

class QuadEdgeLinearNodes : public LinearEdgeReturn
{
protected:
	static constexpr unsigned _edge_nodes[4][2] = {
		{0, 1}, // Edge 0
		{1, 2}, // Edge 1
		{2, 3}, // Edge 2
		{0, 3}  // Edge 3
	};
};

class QuadEdgeQuadraticNodes: public QuadraticEdgeReturn
{
protected:
	static constexpr unsigned _edge_nodes[4][3] = {
		{0, 1, 4}, // Edge 0
		{1, 2, 5}, // Edge 1
		{2, 3, 6}, // Edge 2
		{0, 3, 7}  // Edge 3
	};
};

class TriEdgeLinearNodes: public LinearEdgeReturn
{
protected:
	static constexpr unsigned _edge_nodes[3][2] = {
		{0, 1}, // Edge 0
		{1, 2}, // Edge 1
		{0, 2}  // Edge 2
	};
};

class TriEdgeQuadraticNodes: public QuadraticEdgeReturn
{
protected:
	static constexpr unsigned _edge_nodes[3][3] = {
		{0, 1, 3}, // Edge 0
		{1, 2, 4}, // Edge 1
		{0, 2, 5}, // Edge 2
	};
};

} // end detail

} // end MeshLib


#endif /* EDGE_RULES_H_ */

