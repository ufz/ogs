/**
 * Element.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <vector>
#include "MshEnums.h"

namespace MeshLib {

class Node;

/**
 * Virtual base class for mesh elements.
 */
class Element
{
public:
	const Node* getNode(size_t i) const;
	Node* const* getNodes() const { return _nodes; };

	virtual size_t getDimension() const = 0;

	virtual size_t getNNodes() const = 0;

	size_t getNodeIndex(size_t i) const;

	MshElemType::type getType() const { return _type; };

	size_t getValue() const { return _value; };

	virtual ~Element();

protected:
	Element(Node** nodes, MshElemType::type type, size_t value = 0);
	Element(MshElemType::type type, size_t value = 0);

	MshElemType::type _type;
	size_t _value;
	Node** _nodes;
	std::vector<Element*> _neighbors;

private:

}; /* class */

} /* namespace */

#endif /* ELEMENT_H_ */
