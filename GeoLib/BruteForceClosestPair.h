/**
 * \file BruteForceClosestPair.h
 *
 * Created on 2011-01-25 by Thomas Fischer
 */

#ifndef BRUTEFORCECLOSESTPAIR_H_
#define BRUTEFORCECLOSESTPAIR_H_

#include "ClosestPair.h"

namespace GeoLib {

class BruteForceClosestPair : public ClosestPair {
public:
	BruteForceClosestPair(std::vector<GeoLib::Point*> const & pnts, size_t& id0, size_t& id1);
};

} // end namespace GeoLib

#endif /* BRUTEFORCECLOSESTPAIR_H_ */
