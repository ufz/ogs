// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "GMSHFixedMeshDensity.h"

namespace FileIO
{
namespace GMSH
{
GMSHFixedMeshDensity::GMSHFixedMeshDensity(double mesh_density)
    : _mesh_density(mesh_density)
{
}

void GMSHFixedMeshDensity::initialize(
    std::vector<GeoLib::Point const*> const& vec)
{
    // to avoid a warning here:
    (void)(vec);
}

double GMSHFixedMeshDensity::getMeshDensityAtPoint(
    GeoLib::Point const* const /*unused*/) const
{
    return _mesh_density;
}

double GMSHFixedMeshDensity::getMeshDensityAtStation(
    GeoLib::Point const* const /*unused*/) const
{
    return _mesh_density;
}

}  // namespace GMSH
}  // end namespace FileIO
