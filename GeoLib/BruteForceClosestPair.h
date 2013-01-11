/**
 * \file
 * \author Thomas Fischer
 * \date   2011-01-25
 * \brief  Definition of the BruteForceClosestPair class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BRUTEFORCECLOSESTPAIR_H_
#define BRUTEFORCECLOSESTPAIR_H_

#include "ClosestPair.h"

namespace GeoLib {

class BruteForceClosestPair : public ClosestPair {
public:
	BruteForceClosestPair(std::vector<GeoLib::Point*> const & pnts, std::size_t& id0, std::size_t& id1);
};

} // end namespace GeoLib

#endif /* BRUTEFORCECLOSESTPAIR_H_ */
