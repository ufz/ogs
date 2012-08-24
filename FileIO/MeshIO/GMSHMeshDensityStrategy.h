/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file GMSHMeshDensityStrategy.h
 *
 *  Created on 2010-03-05 by Thomas Fischer
 */

#ifndef GMSHMESHDENSITYSTRATEGY_H_
#define GMSHMESHDENSITYSTRATEGY_H_

#include <ostream>
#include <vector>

// GeoLib
#include "Point.h"

namespace FileIO
{
/**
 * virtual base class GMSHMeshDensityStrategy for classes
 * GMSHAdaptiveMeshDensity, GMSHFixedMeshDensity and GMSHNoMeshDensity
 */
class GMSHMeshDensityStrategy
{
public:
	virtual void init(std::vector<GeoLib::Point const*> const&) = 0;
	virtual double getMeshDensityAtPoint(GeoLib::Point const*const) const = 0;
};

} // end namespace


#endif /* GMSHMESHDENSITYSTRATEGY_H_ */
