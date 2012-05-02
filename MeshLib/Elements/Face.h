/**
 * Face.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef FACE_H_
#define FACE_H_

#include "Element.h"

namespace MeshLib {

/**
 * Virtual base class for 2d mesh elements.
 */
class Face : public Element
{
public:
	virtual double getArea() const { return _area; };

	size_t getDimension() const { return 2; };

	virtual ~Face();

protected:
	Face(Node** nodes, MshElemType::type type, size_t value = 0);
	Face(MshElemType::type type, size_t value = 0);

	virtual double calcArea() = 0;

	double _area;

private:
	

}; /* class */

} /* namespace */

#endif /* FACE_H_ */
