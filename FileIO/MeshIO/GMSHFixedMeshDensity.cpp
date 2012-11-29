/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file GMSHFixedMeshDensity.cpp
 *
 *  Created on 2012-03-05 by Thomas Fischer
 */

#include "MeshIO/GMSHFixedMeshDensity.h"

namespace FileIO {

GMSHFixedMeshDensity::GMSHFixedMeshDensity(double mesh_density) :
	_mesh_density(mesh_density)
{
}

void GMSHFixedMeshDensity::init(std::vector<GeoLib::Point const*> const& vec)
{
	// to avoid a warning here:
	(void)(vec);
}

double GMSHFixedMeshDensity::getMeshDensityAtPoint(GeoLib::Point const*const pnt) const
{
	// to avoid a warning here:
	(void)(const_cast<GeoLib::Point const*>(pnt));
	return _mesh_density;
}

} // end namespace FileIO
