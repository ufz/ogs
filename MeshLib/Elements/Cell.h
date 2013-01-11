/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Cell class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CELL_H_
#define CELL_H_

#include "Element.h"
#include "Face.h"

namespace MeshLib {

/**
 * Virtual base class for 3d mesh elements.
 */
class Cell : public Element
{
public:
	/// Returns the length, area or volume of a 1D, 2D or 3D element
	double getContent() const { return _volume; };

	/// Get dimension of the mesh element.
	unsigned getDimension() const { return 3; };

	/// Get the volume of this 3d element.
	virtual double getVolume() const { return _volume; };

	/// Returns true if the cell is somewhere on the mesh surface and false otherwise.
	bool isOnSurface() const;

	/// Destructor
	virtual ~Cell();

	/**
	 * This method is pure virtual and is inherited from class @sa Element.
	 * It has to be implemented in the derived classes of class Cell!
	 * @return a copy of the object
	 */
	virtual Element* clone() const = 0;

protected:
/*
	/// Constructor for a generic mesh element containing an array of mesh nodes.
	Cell(Node** nodes, MshElemType::type type, unsigned value = 0);
*/
	/// Constructor for a generic mesh element without an array of mesh nodes.
	Cell(unsigned value = 0);

	double _volume;

}; /* class */

}

#endif /* CELL_H_ */

