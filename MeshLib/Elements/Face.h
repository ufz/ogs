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

#ifndef FACE_H_
#define FACE_H_

#include "GeoLib/Point.h"

#include "MathLib/Vector3.h"

#include "Element.h"

namespace MeshLib {

/**
 * Virtual base class for 2d mesh elements.
 */
class Face : public Element
{
public:
	/// Get the area of this 2d element.
	virtual double getArea() const { return _content; }

	/// Returns the surface normal of a 2D element.
	MathLib::Vector3 getSurfaceNormal() const;

	/// Destructor
	virtual ~Face();

	/**
	 * Checks if the node order of an element is correct by testing surface normals.
	 * For 2D elements true is returned if the normal points (roughly) upwards.
	 */
	virtual bool testElementNodeOrder() const;

protected:
	/// Constructor
	/// @param value  element value
	/// @param id     element id
	Face(unsigned value, std::size_t id);

private:


}; /* class */

} /* namespace */

#endif /* FACE_H_ */

