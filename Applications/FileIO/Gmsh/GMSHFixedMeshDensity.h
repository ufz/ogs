/**
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GMSHMeshDensityStrategy.h"

namespace FileIO
{
namespace GMSH
{

class GMSHFixedMeshDensity final : public GMSHMeshDensityStrategy
{
public:
    explicit GMSHFixedMeshDensity(double mesh_density);
    void initialize(std::vector<GeoLib::Point const*> const& vec) override;
    double getMeshDensityAtPoint(
        GeoLib::Point const* const /*unused*/) const override;
    double getMeshDensityAtStation(
        GeoLib::Point const* const /*unused*/) const override;
    ~GMSHFixedMeshDensity() override = default;

private:
    double _mesh_density;
};

}  // namespace GMSH
} // end namespace FileIO
