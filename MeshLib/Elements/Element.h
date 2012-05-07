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
#include "Mesh.h"

namespace MeshLib {

class Node;

/**
 * Virtual base class for mesh elements.
 */
class Element
{
public:
	/// Get node with local index i.
	const Node* getNode(size_t i) const;

	/// Get array of element nodes.
	Node* const* getNodes() const { return _nodes; };

	/// Get dimension of the mesh element.
	virtual size_t getDimension() const = 0;

	/// Get the number of nodes for this element.
	virtual size_t getNNodes() const = 0;

	/// Get the global index for the node with local index i.
	size_t getNodeIndex(size_t i) const;

	/// Get the type of the mesh element (as a MshElemType-enum).
	MshElemType::type getType() const { return _type; };

	/// Get the value for this element.
	size_t getValue() const { return _value; };

	/// Destructor
	virtual ~Element();

protected:
	/// Constructor for a generic mesh element containing an array of mesh nodes.
	Element(Node** nodes, MshElemType::type type, size_t value = 0);

	/// Constructor for a generic mesh element without an array of mesh nodes.
	Element(MshElemType::type type, size_t value = 0);

	/**
	 * Get an editale Node.
	 * This method is called by Mesh::addElement(Element*), see friend definition.
	 */
	Node* getNode(size_t i);

	MshElemType::type _type;
	Node** _nodes;
	size_t _value;
	std::vector<Element*> _neighbors;

private:

/* friend functions: */
	friend void Mesh::addElement(Element*);

}; /* class */

} /* namespace */

#endif /* ELEMENT_H_ */

