/*
 * LinkedTriangle.h
 *
 *  Created on: Mar 25, 2010
 *      Author: fischeth
 */

#ifndef LINKEDTRIANGLE_H_
#define LINKEDTRIANGLE_H_

#include <vector>
#include <iostream>

// GEO
#include "Triangle.h"
#include "Point.h"

namespace MathLib {

class LinkedTriangle : public GEOLIB::Triangle {
public:
	LinkedTriangle(std::vector<GEOLIB::Point*> const&pnt_vec, size_t pnt_a,
			size_t pnt_b, size_t pnt_c, LinkedTriangle* tri_a,
			LinkedTriangle* tri_b, LinkedTriangle* tri_c);
	virtual ~LinkedTriangle();

	void setNeighborTriangle (size_t idx, LinkedTriangle* tri);
	void setNeighborTriangleByPointIdx (size_t idx, LinkedTriangle* tri);

	LinkedTriangle* getNeighborTriangle (size_t idx);
	size_t getIdxOfNeighborTriangle(LinkedTriangle* tri);

	size_t getIdxOfPoint (size_t i) const;

	void write (std::ostream &os) const
	{
		os << _pnt_ids[0] << " " << _pnt_ids[1] << " " << _pnt_ids[2];
	}

	void writeNeighbor (std::ostream &os, size_t idx) const;

private:
	/**
	 * pointers to the neighbor triangles
	 */
	LinkedTriangle* _neighbor_triangles[3];
};

/** overload the output operator for class LinkedTriangle */
std::ostream& operator<< (std::ostream &os, const LinkedTriangle &tri);

} // end namespace MathLib

#endif /* LINKEDTRIANGLE_H_ */
