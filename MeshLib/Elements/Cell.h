/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Cell class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CELL_H_
#define CELL_H_

#include "Element.h"

namespace MeshLib {

/**
 * Virtual base class for 3d mesh elements.
 */
class Cell : public Element
{
public:
	/// Get the volume of this 3d element.
	virtual double getVolume() const { return getContent(); }

	/// Destructor
	virtual ~Cell();

	/**
	 * Checks if the node order of an element is correct by testing surface normals.
	 * For 3D elements true is returned if the normals of all faces points away from the centre of 
	 * the element.
	 * Note: This method might give wrong results if something else is wrong with the element 
	 * (non-planar faces, non-convex geometry, possibly zero volume) which causes the calculated 
	 * center of gravity to lie outside of the actual element
	 */
	virtual bool testElementNodeOrder() const;

protected:
/*
	/// Constructor for a generic mesh element containing an array of mesh nodes.
	Cell(Node** nodes, MeshElemType type, unsigned value = 0);
*/
	/// Constructor for a generic mesh element without an array of mesh nodes.
	Cell(unsigned value = 0, std::size_t id = std::numeric_limits<std::size_t>::max());
}; /* class */

}

#endif /* CELL_H_ */

