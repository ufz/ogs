/*
 * BruteForceClosestPair.h
 *
 *  Created on: Jan 25, 2011
 *      Author: TF
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
