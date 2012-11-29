/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file GMSHFixedMeshDensity.h
 *
 *  Created on 2012-03-05 by Thomas Fischer
 */

#ifndef GMSHFIXEDMESHDENSITY_H_
#define GMSHFIXEDMESHDENSITY_H_

#include "GMSHMeshDensityStrategy.h"

namespace FileIO {

class GMSHFixedMeshDensity : public GMSHMeshDensityStrategy
{
public:
	GMSHFixedMeshDensity(double mesh_density);
	void init(std::vector<GeoLib::Point const*> const& vec);
	double getMeshDensityAtPoint(GeoLib::Point const*const) const;

private:
	double _mesh_density;
};

}

#endif /* GMSHFIXEDMESHDENSITY_H_ */
