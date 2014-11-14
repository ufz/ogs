/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Face class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef FACE_H_
#define FACE_H_

#include <limits>

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
	/// Constant: Dimension of this mesh element
	static const unsigned dimension;

	/// Get the area of this 2d element.
	virtual double getArea() const { return _area; };

	/// Returns the length, area or volume of a 1D, 2D or 3D element
	double getContent() const { return _area; };

	/// Get dimension of the mesh element.
	unsigned getDimension() const { return dimension; };

	/// Returns the face i of the element.
	const Element* getFace(unsigned i) const { return this->getEdge(i); };

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned /*i*/) const { return 2; };

	/// 2D elements have no faces.
	unsigned getNFaces() const { return 0; };

	/// Returns the surface normal of a 2D element.
	MathLib::Vector3 getSurfaceNormal() const;

	/// Destructor
	virtual ~Face();

	/**
	 * This method is pure virtual and is inherited from class @sa Element.
	 * It has to be implemented in the derived classes of class Face!
	 * @return a copy of the object
	 */
	virtual Element* clone() const = 0;

	/**
	 * Checks if the node order of an element is correct by testing surface normals.
	 * For 2D elements true is returned if the normal points (roughly) upwards.
	 */
	virtual bool testElementNodeOrder() const;

protected:
/*
	/// Constructor for a generic mesh element containing an array of mesh nodes.
	Face(Node** nodes, MeshElemType type, unsigned value = 0);
*/
	/// Constructor for a generic mesh element without an array of mesh nodes.
	Face(unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());

	double _area;

private:


}; /* class */

} /* namespace */

#endif /* FACE_H_ */

