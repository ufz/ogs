/**
 * FemElem.h
 *
 *      Date: 2012/05/03
 *      Author: KR
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
	const GEOLIB::Point& getCentreOfGravity() const { return _centre_of_gravity; };

protected:
	/// Constructor.
	FemElem();

	/// Copy constructor
	FemElem(const FemElem &elem);

	/// Calculate centre of gravity
	virtual void calcCoG() = 0;

	GEOLIB::Point _centre_of_gravity;

}; /* class */

} /* namespace */

#endif /* FEMELEM_H_ */

