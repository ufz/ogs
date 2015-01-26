/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Face class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EDGE_H_
#define EDGE_H_

#include <limits>

#include "Element.h"

namespace MeshLib {

/**
 * Virtual base class for 1d mesh elements.
 */
class Edge : public Element
{
public:
	/// Get the length of this 1d element.
	virtual double getLength() const { return _content; }

	/// Destructor
	virtual ~Edge() {}

	/**
		* Checks if the node order of an element is correct by testing surface normals.
		* For 1D elements this always returns true.
		*/
	virtual bool testElementNodeOrder() const { return true; }

protected:
	/// Constructor for a generic mesh element without an array of mesh nodes.
	Edge(unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max()) : Element(value, id) {}


}; /* class */

} /* namespace */

#endif /* EDGE_H_ */

