/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Tri.h
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#ifndef TRI_H_
#define TRI_H_

#include "Edge.h"
#include "Node.h"
#include "Face.h"

#include "MathTools.h"


namespace MeshLib {

/**
 * This class represents a 2d triangle element. The following sketch shows the node and edge numbering.
 * @anchor TriNodeAndEdgeNumbering
 * @code
 *
 *          2
 *          o
 *         / \
 *        /   \
 *      2/     \1
 *      /       \
 *     /         \
 *    0-----------1
 *          0
 *
 * @endcode
 */
template <unsigned ORDER, unsigned NNODES>
class TemplateTri : public Face
{
public:
	/// Constructor with an array of mesh nodes.
	TemplateTri(Node* nodes[NNODES], unsigned value = 0) : Face(value)
	{
		_nodes = nodes;
		_neighbors = new Element*[3];
		for (unsigned i=0; i<3; i++)
			_neighbors[i] = NULL;
		this->_area = this->computeVolume();
	}

	/// Copy constructor
	TemplateTri(const TemplateTri<ORDER, NNODES> &tri) : Face(tri.getValue())
	{
		_nodes = new Node*[NNODES];
		for (unsigned i=0; i<NNODES; i++)
		{
			_nodes[i] = tri._nodes[i];
		}

		_neighbors = new Element*[3];
		for (unsigned i=0; i<3; i++) {
			_neighbors[i] = tri._neighbors[i];
		}

		_area = tri.getArea();
	}

	/// Destructor
	virtual ~TemplateTri() {};

	/// Get the number of edges for this element.
	unsigned getNEdges() const { return 3; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 3; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes(unsigned order = 1) const
	{
		return order == ORDER ? NNODES : 3;
	}

	/**
	 * Method returns the type of the element. In this case TRIANGLE will be returned.
	 * @return MshElemType::TRIANGLE
	 */
	virtual MshElemType::type getType() const { return MshElemType::TRIANGLE; }

	/// Returns true if these two indices form an edge and false otherwise
	bool isEdge(unsigned idx1, unsigned idx2) const
	{
		for (unsigned i(0); i<3; i++)
		{
			if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
			if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
		}
		return false;
	}

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the TemplateTri instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const
	{
		return new TemplateTri<ORDER, NNODES>(*this);
	}


	/**
	 * This method should be called after at least two nodes of the triangle
	 * element are collapsed. As a consequence of the node collapsing an edge
	 * of the triangle will be collapsed. If one of the edges is collapsed we
	 * obtain an edge. In this case the method will create the appropriate
	 * object of class Edge.
	 * @return an Edge object or NULL
	 */
	virtual Element* reviseElement() const
	{
		// try to create an edge
		if (_nodes[0] == _nodes[1] || _nodes[1] == _nodes[2]) {
			Node** nodes (new Node*[2]);
			nodes[0] = _nodes[0];
			nodes[1] = _nodes[2];
			return new Edge(nodes, _value);
		}

		if (_nodes[0] == _nodes[2]) {
			Node** nodes (new Node*[2]);
			nodes[0] = _nodes[0];
			nodes[1] = _nodes[1];
			return new Edge(nodes, _value);
		}

		return NULL;
	}

protected:
	/// Calculates the area of the triangle by returning half of the area of the corresponding parallelogram.
	double computeVolume()
	{
		return MathLib::calcTriangleArea(_nodes[0]->getCoords(), _nodes[1]->getCoords(), _nodes[2]->getCoords());
	}

protected:
	/// Return a specific edge node.
	inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { return _nodes[_edge_nodes[edge_id][node_id]]; };

	/// Returns the ID of a face given an array of nodes.
	unsigned identifyFace(Node* nodes[3]) const
	{
		for (unsigned i=0; i<3; i++)
		{
			unsigned flag(0);
			for (unsigned j=0; j<2; j++)
				for (unsigned k=0; k<2; k++)
					if (_nodes[_edge_nodes[i][j]] == nodes[k])
						flag++;
			if (flag==2)
				return i;
		}
		return std::numeric_limits<unsigned>::max();
	}

	static const unsigned _edge_nodes[3][2];
}; /* class */

typedef TemplateTri<1,3> Tri;

template <unsigned ORDER, unsigned NNODES>
const unsigned TemplateTri<ORDER, NNODES>::_edge_nodes[3][2] = {
		{0, 1}, // Edge 0
		{1, 2}, // Edge 1
		{0, 2}  // Edge 2
	};

} /* namespace */

#endif /* TRI_H_ */

