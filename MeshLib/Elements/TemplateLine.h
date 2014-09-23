/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Line class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEMPLATELINE_H_
#define TEMPLATELINE_H_

#include <array>
#include <limits>
#include <cmath>

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
template<unsigned NNODES, CellType CELLLINETYPE>
class TemplateLine : public Element
{
public:
	/// Constant: Dimension of this mesh element
	static const unsigned dimension = 1u;

	/// Constant: The number of all nodes for this element
	static const unsigned n_all_nodes = NNODES;

	/// Constant: The number of base nodes for this element
	static const unsigned n_base_nodes = 2u;

	/// Constructor with an array of mesh nodes.
	TemplateLine(Node* nodes[NNODES], unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Constructs a line from array of Node pointers.
	TemplateLine(std::array<Node*, NNODES> const& nodes, unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	/// Copy constructor
	TemplateLine(const TemplateLine &line);

	/// Destructor
	virtual ~TemplateLine();

	/// Compute the minimum and maximum squared edge length for this element
	void computeSqrEdgeLengthRange(double &min, double &max) const { min = _length; max = _length; };

	/// Returns the length, area or volume of a 1D, 2D or 3D element
	double getContent() const { return _length; };

	/// Get dimension of the mesh element.
	unsigned getDimension() const { return dimension; };

	/// Returns the edge i of the element.
	const Element* getEdge(unsigned /*i*/) const { return nullptr; };

	/// Returns the face i of the element.
	const Element* getFace(unsigned /*i*/) const { return nullptr; };

	/// Get the length of this 1d element.
	double getLength() const { return _length; };

	/// 1D elements have no edges
	unsigned getNEdges() const { return 1; };

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned /*i*/) const { return 0; };

	/// Get the number of faces for this element.
	unsigned getNFaces() const { return 0; };

	/// Get the number of neighbors for this element.
	unsigned getNNeighbors() const { return 2; }

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes(bool all = false) const
	{
		return all ? n_all_nodes : n_base_nodes;
	}

	/**
		* Method returns the type of the element. In this case LINE will be returned.
		* @return MeshElemType::LINE
		*/
	virtual MeshElemType getGeomType() const { return MeshElemType::LINE; }

	/**
		* Get the type of the element in context of the finite element method.
		* @return a value of the enum FEMElemType::type
		*/
	virtual CellType getCellType() const { return CELLLINETYPE; }

	/// Returns true if these two indices form an edge and false otherwise
	bool isEdge(unsigned idx1, unsigned idx2) const
	{
		if (0==idx1 && 1==idx2) return true;
		if (1==idx1 && 0==idx2) return true;
		return false;
	}

    /// Returns true if pnt is located on the line segment and false otherwise
    bool isPntInElement(GeoLib::Point const& pnt) const;

	/**
		* Tests if the element is geometrically valid.
		* @param check_zero_volume indicates if area == 0 should be checked
		* @return error code (0 = okay, 1 = zero volume)
		*/
	virtual ElementErrorCode validate() const;

	/**
		* Checks if the node order of an element is correct by testing surface normals.
		* For 1D elements this always returns true.
		*/
	virtual bool testElementNodeOrder() const { return true; }

	/**
		* Method clone is inherited from class Element. It makes a deep copy of the TemplateLine instance.
		* @return an exact copy of the object
		*/
	virtual Element* clone() const
	{
		return new TemplateLine<NNODES,CELLLINETYPE>(*this);
	}

protected:
	double computeVolume()
	{
		return sqrt(MathLib::sqrDist(_nodes[0]->getCoords(), _nodes[1]->getCoords()));
	}

	/// Returns the specified node.
	Node* getEdgeNode(unsigned edge_id, unsigned node_id) const 
	{ 
		if (edge_id==0 && node_id<2)
			return _nodes[node_id];
		return nullptr;
	};

	/// Returns the ID of a face given an array of nodes.
	/// As element faces are always element->getDimensionality() - 1, the "face" of a line is just a node
	/// and the method returns if another line is connected to node[0] or node[1] of the line.
	unsigned identifyFace(Node* nodes[1]) const
	{
		if (nodes[0] == _nodes[0])
			return 0;
		if (nodes[0] == _nodes[1])
			return 1;
		return std::numeric_limits<unsigned>::max();
	}

	double _length;

}; /* class */

} /* namespace */

#include "TemplateLine-impl.h"

#endif /* TEMPLATELINE_H_ */

