/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
 *
 *
 * \file Cell.h
 *
 * Created on 2012-05-02 by Karsten Rink
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
	/// Returns the length, area or volume of a 1D, 2D or 3D element
	double getContent() const { return _volume; };

	/// Get dimension of the mesh element.
	unsigned getDimension() const { return 3; };

	/// Get the volume of this 3d element.
	virtual double getVolume() const { return _volume; };

	/// Destructor
	virtual ~Cell();

protected:
/*
	/// Constructor for a generic mesh element containing an array of mesh nodes.
	Cell(Node** nodes, MshElemType::type type, unsigned value = 0);
*/
	/// Constructor for a generic mesh element without an array of mesh nodes.
	Cell(MshElemType::type type, unsigned value = 0);

	/// Calculate the volume of this 3d element.
	virtual double computeVolume() = 0;

	double _volume;

}; /* class */

}

#endif /* CELL_H_ */

