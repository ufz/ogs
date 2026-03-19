// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

#include "DiameterProfile.h"

namespace BaseLib
{
class ConfigTree;
}
namespace MeshLib
{
class Node;
}
namespace ParameterLib
{
struct ParameterBase;
}
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct BoreholeGeometry
{
    /**
     * Total length/depth of the BHE
     * unit is m
     */
    double const length;

    /**
     * Z-coordinate of the wellhead (topmost BHE node) [m].
     * Used to compute distance from wellhead for multi-section lookup.
     */
    double const wellhead_z;

    /// Section boundaries and borehole diameters.
    DiameterProfile const sections;
};

BoreholeGeometry createBoreholeGeometry(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::vector<MeshLib::Node*> const& bhe_nodes);

}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
