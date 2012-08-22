/*
 * GMSHMeshDensityStrategy.h
 *
 *  Created on: Mar 5, 2012
 *      Author: TF
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
