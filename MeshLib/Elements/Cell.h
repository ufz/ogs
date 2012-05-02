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
	virtual double getVolume() const { return _volume; };

	size_t getDimension() const { return 3; };

	virtual ~Cell();

protected:
	Cell(Node** nodes, MshElemType::type type, size_t value = 0);
	Cell(MshElemType::type type, size_t value = 0);

	virtual double calcVolume() = 0;

	double _volume;

}; /* class */

} /* namespace */

#endif /* CELL_H_ */
