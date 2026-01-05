// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

namespace GeoLib
{
class Point;
}

namespace FileIO
{
namespace GMSH
{
/**
 * virtual base class GMSHMeshDensityStrategy for classes
 * GMSHAdaptiveMeshDensity and GMSHFixedMeshDensity.
 */
class GMSHMeshDensityStrategy
{
public:
    virtual ~GMSHMeshDensityStrategy() = default;
    virtual void initialize(std::vector<GeoLib::Point const*> const&) = 0;
    virtual double getMeshDensityAtPoint(GeoLib::Point const*const) const = 0;
    virtual double getMeshDensityAtStation(GeoLib::Point const*const) const = 0;
};

}  // end namespace GMSH
}  // end namespace FileIO
