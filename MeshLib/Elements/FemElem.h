/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-03
 * \brief  Definition of the FemElem class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef FEMELEM_H_
#define FEMELEM_H_

#include "Point.h"

namespace MeshLib {

/**
 * Virtual base class for elements of finite element meshes.
 */
class FemElem
{
public:
	/// Destructor
	virtual ~FemElem();

	/// Get the number of nodes for this element.
	const GeoLib::Point& getCentroid() const { return _centroid; };

protected:
	/// Constructor.
	FemElem();

	/// Copy constructor
	FemElem(const FemElem &elem);

	/// Calculate centre of gravity
	virtual void calcCentroid() = 0;

	GeoLib::Point _centroid;

}; /* class */

} /* namespace */

#endif /* FEMELEM_H_ */

