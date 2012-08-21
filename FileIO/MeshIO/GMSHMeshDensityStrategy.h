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

// GEOLIB
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
	virtual void init(std::vector<GEOLIB::Point const*> const&) = 0;
	virtual double getMeshDensityAtPoint(GEOLIB::Point const*const) const = 0;
};

} // end namespace


#endif /* GMSHMESHDENSITYSTRATEGY_H_ */
