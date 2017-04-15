/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
