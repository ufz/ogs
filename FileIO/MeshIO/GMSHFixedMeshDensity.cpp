/*
 * GMSHFixedMeshDensity.cpp
 *
 *  Created on: Mar 5, 2012
 *      Author: TF
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
