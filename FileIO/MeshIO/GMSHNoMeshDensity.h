/*
 * GMSHNoMeshDensity.h
 *
 *  Created on: Mar 5, 2012
 *      Author: fischeth
 */

#ifndef GMSHNOMESHDENSITY_H_
#define GMSHNOMESHDENSITY_H_

#include "GMSHMeshDensityStrategy.h"

namespace FileIO {

class GMSHNoMeshDensity: public FileIO::GMSHMeshDensityStrategy {
public:
	GMSHNoMeshDensity() {};
	void init(std::vector<GeoLib::Point const*> const& vec)
	{
		// to avoid a warning here:
		(void)(vec);
	}

	double getMeshDensityAtPoint(GeoLib::Point const*const pnt) const
	{
		// to avoid a warning here:
		(void)(pnt);
		return 0.0;
	}
};

}

#endif /* GMSHNOMESHDENSITY_H_ */
