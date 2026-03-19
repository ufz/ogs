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
template <typename T>
struct Parameter;
}  // namespace ParameterLib
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

    /// Section boundaries and borehole diameters.
    DiameterProfile const sections;

    /// Stored for rebuilding geometry with different BHE nodes
    /// (grouped BHE definitions).  May be nullptr for test-constructed objects.
    ParameterLib::Parameter<double> const* diameter_param{nullptr};

    /// Rebuild this geometry for different BHE nodes.
    /// Used for grouped BHE definitions where config is shared but each BHE
    /// has its own node set.
    BoreholeGeometry rebuildForNodes(
        std::vector<MeshLib::Node*> const& new_nodes) const;
};

BoreholeGeometry createBoreholeGeometry(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::vector<MeshLib::Node*> const& bhe_nodes);

}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
