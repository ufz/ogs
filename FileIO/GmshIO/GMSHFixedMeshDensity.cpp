/**
 * \file
 * \author Thomas Fischer
 * \date   2012-03-05
 * \brief  Implementation of the GMSHFixedMeshDensity class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GmshIO/GMSHFixedMeshDensity.h"

namespace FileIO 
{
namespace GMSH {

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

}
} // end namespace FileIO
