/*
 * ClosestPair.h
 *
 *  Created on: Jan 25, 2011
 *      Author: TF
 */

#ifndef CLOSESTPAIR_H_
#define CLOSESTPAIR_H_

// STL
#include <vector>

// GEOLIB
#include "Point.h"

namespace GEOLIB {

class ClosestPair
{
public:
	ClosestPair (std::vector<GEOLIB::Point*> const & pnts, size_t id0, size_t id1) :
		_pnts (pnts), _id0 (id0), _id1 (id1)
	{}

protected:
	std::vector<GEOLIB::Point*> const & _pnts;
	size_t _id0;
	size_t _id1;
};

} // end namespace GEOLIB

#endif /* CLOSESTPAIR_H_ */
