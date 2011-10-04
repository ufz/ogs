/*
 * LinkedTriangle.cpp
 *
 *  Created on: Mar 25, 2010
 *      Author: fischeth
 */

#include "LinkedTriangle.h"

namespace MathLib {

LinkedTriangle::LinkedTriangle(std::vector<GEOLIB::Point*> const &pnt_vec,
		size_t pnt_a, size_t pnt_b, size_t pnt_c, LinkedTriangle* tri_a,
		LinkedTriangle* tri_b, LinkedTriangle* tri_c) :
	GEOLIB::Triangle(pnt_vec, pnt_a, pnt_b, pnt_c)
{
	if (tri_a) _neighbor_triangles[0] = tri_a;
	else _neighbor_triangles[0] = NULL;
	if (tri_b) _neighbor_triangles[1] = tri_b;
	else _neighbor_triangles[1] = NULL;
	if (tri_c) _neighbor_triangles[2] = tri_c;
	else _neighbor_triangles[2] = NULL;
}

void LinkedTriangle::setNeighborTriangle(size_t idx, LinkedTriangle* tri)
{
	assert(idx < 3);
	_neighbor_triangles[idx] = tri;
}

void LinkedTriangle::setNeighborTriangleByPointIdx (size_t idx, LinkedTriangle* tri)
{
	if (idx == _pnt_ids[0]) setNeighborTriangle (0, tri);
	else {
		if (idx == _pnt_ids[1]) setNeighborTriangle (1, tri);
		else setNeighborTriangle (2, tri);
	}
}

LinkedTriangle* LinkedTriangle::getNeighborTriangle(size_t idx)
{
	assert(idx < 3);
	return _neighbor_triangles[idx];
}

size_t LinkedTriangle::getIdxOfNeighborTriangle(LinkedTriangle* tri)
{
	 if (tri == _neighbor_triangles[0]) return 0;
	 if (tri == _neighbor_triangles[1]) return 1;
	 if (tri == _neighbor_triangles[2]) return 2;
	 return 3;
}

size_t LinkedTriangle::getIdxOfPoint (size_t i) const
{
	if (_pnt_ids[0] == i) return 0;
	if (_pnt_ids[1] == i) return 1;
	if (_pnt_ids[2] == i) return 2;
	return 3;
}

void LinkedTriangle::writeNeighbor (std::ostream &os, size_t idx) const
{
	if (_neighbor_triangles[idx]) {
		os << "[";
		(_neighbor_triangles[idx])->write (os);
		os << "]";
	}
	else os << "NULL";
}

LinkedTriangle::~LinkedTriangle()
{
	// TODO Auto-generated destructor stub
}

std::ostream& operator<< (std::ostream &os, const LinkedTriangle &tri)
{
	tri.write (os);
	return os;
}

} // end namespace MathLib
