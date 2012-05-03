/**
 * Node.h
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#ifndef NODE_H_
#define NODE_H_

#include <cstdlib>
#include <limits>
#include <vector>

#include "PointWithID.h"
#include "Mesh.h"

namespace MeshLib {

class Element;

/**
 * A mesh node with coordinates in 3D space.
 */
class Node : public GEOLIB::PointWithID
{
public:
	/// Constructor using a coordinate array
	Node(const double coords[3], size_t id = std::numeric_limits<size_t>::max());

	/// Constructor using single coordinates
	Node(double x, double y, double z, size_t id = std::numeric_limits<size_t>::max());

	/// Copy constructor
	Node(const Node &node);

	/// Get an element the node is part of.
	const Element* getElement(size_t idx) const { return _elements[idx]; };

	/// Get all elements the node is part of.
	const std::vector<const Element*> getElements() const { return _elements; };

	/// Get number of elements the node is part of.
	size_t getNElements() const { return _elements.size(); };

	/// Destructor
	~Node();

protected:
	/**
	 * Add an element the node is part of.
	 * This method is called by Mesh::addElement(Element*), see friend definition.
	 */
	void addElement(const Element* elem) { _elements.push_back(elem); };

	/// Remove an element the node is part of.
	//void removeElement(const Element* elem);

	std::vector<const Element*> _elements;

/* friend functions: */
	friend void Mesh::addElement(Element*);
	//friend void Mesh::removeElement(size_t);

}; /* class */

} /* namespace */

#endif /* NODE_H_ */

