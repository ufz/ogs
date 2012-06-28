/**
 * \file ClosestPair.h
 *
 * Created on 2011-01-25 by Thomas Fischer
 */

#ifndef CLOSESTPAIR_H_
#define CLOSESTPAIR_H_

// STL
#include <vector>

// GeoLib
#include "Point.h"

namespace GeoLib {

class ClosestPair
{
public:
	ClosestPair (std::vector<GeoLib::Point*> const & pnts, size_t id0, size_t id1) :
		_pnts (pnts), _id0 (id0), _id1 (id1)
	{}

protected:
	std::vector<GeoLib::Point*> const & _pnts;
	size_t _id0;
	size_t _id1;
};

} // end namespace GeoLib

#endif /* CLOSESTPAIR_H_ */
