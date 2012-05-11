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
	/// Get the area of this 2d element.
	virtual double getArea() const { return _area; };

	/// Get dimension of the mesh element.
	unsigned getDimension() const { return 2; };

	/// Returns the face i of the element.
	const Element* getFace(unsigned i) const { return this->getEdge(i); };

	/// Get the number of nodes for face i.
	unsigned getNFaceNodes(unsigned i) const { return 2; };

	/// 2D elements have no faces.
	unsigned getNFaces() const { return 0; };

	/// Destructor
	virtual ~Face();

protected:
/*
	/// Constructor for a generic mesh element containing an array of mesh nodes.
	Face(Node** nodes, MshElemType::type type, unsigned value = 0);
*/
	/// Constructor for a generic mesh element without an array of mesh nodes.
	Face(MshElemType::type type, unsigned value = 0);

	/// Calculate the area of this 2d element.
	virtual double computeArea() = 0;

	double _area;

private:
	

}; /* class */

} /* namespace */

#endif /* FACE_H_ */

