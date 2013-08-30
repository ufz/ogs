/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Edge class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATEEDGE_H_
#define TEMPLATEEDGE_H_

#include <array>
#include <limits>

#include "MeshEnums.h"
#include "Element.h"
#include "Node.h"

#include "MathTools.h"


namespace MeshLib {

/**
 * A 1d Edge or Line Element.
 * @code
 *  0--------1
 * @endcode
 */
template<unsigned NNODES, CellType CELLEDGETYPE>
class TemplateEdge : public Element
{
public:
	/// Constructor with an array of mesh nodes.
	TemplateEdge(Node* nodes[NNODES], unsigned value = 0);

	/// Constructs an edge from array of Node pointers.
	TemplateEdge(std::array<Node*, NNODES> const& nodes, unsigned value = 0);

	/// Copy constructor
	TemplateEdge(const TemplateEdge &edge);

	/// Destructor
	virtual ~TemplateEdge();

	/// Returns the length, area or volume of a 1D, 2D or 3D element
	double getContent() const { return _length; };

	/// Returns the edge i of the element.
	const Element* getEdge(unsigned i) const { (void)i; return NULL; };

	/// Returns the face i of the element.
	const Element* getFace(unsigned i) const { (void)i; return NULL; };

	/// Compute the minimum and maximum squared edge length for this element
	void computeSqrEdgeLengthRange(double &min, double &max) const { min = _length; max = _length; };

	/// 1D elements have no edges
	unsigned getNEdges() const { return 0; };

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned i) const { (void)i; return 0; };

	/// Get the number of faces for this element.
	unsigned getNFaces() const { return 0; };

	/// Get the length of this 1d element.
	double getLength() const { return _length; };

	/// Get dimension of the mesh element.
	unsigned getDimension() const { return 1; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 0; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes(bool all = false) const
	{
		return all ? NNODES : 2;
	}

	/**
	 * Method returns the type of the element. In this case EDGE will be returned.
	 * @return MeshElemType::EDGE
	 */
	virtual MeshElemType getGeomType() const { return MeshElemType::EDGE; }

	/**
	 * Get the type of the element in context of the finite element method.
	 * @return a value of the enum FEMElemType::type
	 */
	virtual CellType getCellType() const { return CELLEDGETYPE; }

	/// Returns true if these two indices form an edge and false otherwise
	bool isEdge(unsigned idx1, unsigned idx2) const
	{
		if (0==idx1 && 1==idx2) return true;
		if (1==idx1 && 0==idx2) return true;
		return false;
	}

	virtual Element* clone() const
	{
		return new TemplateEdge<NNODES,CELLEDGETYPE>(*this);
	}

	/**
	 * If for instance two nodes of the element are collapsed the Edge element disappears.
	 * @return NULL
	 */
	virtual Element* reviseElement() const
	{
		if (_nodes[0] == _nodes[1]) {
			return NULL;
		}

		return NULL;
	}

protected:
	double computeVolume()
	{
		return sqrt(MathLib::sqrDist(_nodes[0]->getCoords(), _nodes[1]->getCoords()));
	}

	/// 1D elements have no edges.
	Node* getEdgeNode(unsigned edge_id, unsigned node_id) const { (void)edge_id; (void)node_id; return NULL; };

	/// 1D elements have no faces.
	Node* getFaceNode(unsigned face_id, unsigned node_id) const { (void)face_id; (void)node_id; return NULL; };

	/// Returns the ID of a face given an array of nodes (but is not applicable for edges!).
	unsigned identifyFace(Node* [3]/*nodes[3]*/) const { return std::numeric_limits<unsigned>::max(); };

	double _length;

}; /* class */

} /* namespace */

#include "TemplateEdge.tpp"

#endif /* TEMPLATEEDGE_H_ */

