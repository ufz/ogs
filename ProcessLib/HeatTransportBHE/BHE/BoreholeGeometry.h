// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
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

    /// Borehole diameter parameter [m].
    /// May be a ConstantParameter or a MeshElementParameter for
    /// spatially varying diameter.  Evaluated on-the-fly via SpatialPosition.
    ParameterLib::Parameter<double> const& diameter;
};

BoreholeGeometry createBoreholeGeometry(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters);

}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
