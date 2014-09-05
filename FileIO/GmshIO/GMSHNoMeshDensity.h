/**
 * \file
 * \author Thomas Fischer
 * \date   Mar 5, 2012
 * \brief  Definition of the GMSHNoMeshDensity class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GMSHNOMESHDENSITY_H_
#define GMSHNOMESHDENSITY_H_

#include "GMSHMeshDensityStrategy.h"

namespace FileIO 
{
namespace GMSH {

class GMSHNoMeshDensity: public FileIO::GMSH::GMSHMeshDensityStrategy {
public:
	GMSHNoMeshDensity() {};
	virtual ~GMSHNoMeshDensity() {};
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
}

#endif /* GMSHNOMESHDENSITY_H_ */
