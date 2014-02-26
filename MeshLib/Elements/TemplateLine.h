/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Line class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "Edge.h"
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
class TemplateLine : public Edge
{
public:
	/// Constructor with an array of mesh nodes.
	TemplateLine(Node* nodes[NNODES], unsigned value = 0);

	/// Constructs a line from array of Node pointers.
	TemplateLine(std::array<Node*, NNODES> const& nodes, unsigned value = 0);

	/// Copy constructor
	TemplateLine(const TemplateLine &line);

	/// Destructor
	virtual ~TemplateLine();

	/// Compute the minimum and maximum squared edge length for this element
	void computeSqrEdgeLengthRange(double &min, double &max) const { min = _length; max = _length; };

	/// Get the number of nodes for this element.
	virtual unsigned getNNodes(bool all = false) const
	{
		return all ? NNODES : 2;
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

	/**
	 * Tests if the element is geometrically valid.
	 * @param check_zero_volume indicates if area == 0 should be checked
	 * @return error code (0 = okay, 1 = zero volume)
	 */
	virtual ElementErrorCode validate() const;

	/**
	 * Method clone is inherited from class Element. It makes a deep copy of the TemplateLine instance.
	 * @return an exact copy of the object
	 */
	virtual Element* clone() const
	{
		return new TemplateLine<NNODES,CELLLINETYPE>(*this);
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

}; /* class */

} /* namespace */

#include "TemplateLine.tpp"

#endif /* TEMPLATELINE_H_ */

