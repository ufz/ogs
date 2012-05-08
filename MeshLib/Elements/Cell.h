/**
 * Cell.h
 *
 *      Date: 2012/05/02
 *      Author: KR
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
	virtual double getVolume() const { return _volume; };

	/// Get dimension of the mesh element.
	unsigned getDimension() const { return 3; };

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
	virtual double calcVolume() = 0;

	double _volume;

}; /* class */

}		

#endif /* CELL_H_ */

