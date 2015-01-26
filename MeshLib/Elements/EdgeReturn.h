/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EDGERETURN_H_
#define EDGERETURN_H_

#include "logog/include/logog.hpp"

namespace MeshLib
{

class Element;

/// Returns always null pointer
class DummyEdgeReturn
{
public:
	/// Returns i-th edge of the given element
	static const Element* getEdge(const Element* /*e*/, unsigned /*i*/)
	{
		return nullptr;
	}
};

/// Returns linear order edge
class LinearEdgeReturn
{
public:
	/// Returns i-th edge of the given element
	static const Element* getEdge(const Element* e, unsigned i);
};

/// Returns quadratic order edge
class QuadraticEdgeReturn
{
public:
	/// Returns i-th edge of the given element
	static const Element* getEdge(const Element* e, unsigned i);
};

} // end MeshLib


#endif /* EDGERETURN_H_ */

