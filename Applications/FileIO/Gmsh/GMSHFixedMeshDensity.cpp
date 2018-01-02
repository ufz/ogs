/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GMSHFixedMeshDensity.h"

namespace FileIO
{
namespace GMSH
{

GMSHFixedMeshDensity::GMSHFixedMeshDensity(double mesh_density) :
    _mesh_density(mesh_density)
{
}

void GMSHFixedMeshDensity::initialize(std::vector<GeoLib::Point const*> const& vec)
{
    // to avoid a warning here:
    (void)(vec);
}

double GMSHFixedMeshDensity::getMeshDensityAtPoint(GeoLib::Point const*const) const
{
    return _mesh_density;
}

double GMSHFixedMeshDensity::getMeshDensityAtStation(GeoLib::Point const*const) const
{
    return _mesh_density;
}

}
} // end namespace FileIO
