/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www./**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
opengeosys.com/LICENSE.txt
 *
 *
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
