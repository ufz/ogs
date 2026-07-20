// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <numbers>
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
struct Pipe
{
    /// Inner diameter [m].
    double const diameter;

    /// Wall thickness [m].
    double const wall_thickness;

    /// Wall thermal conductivity [W/(m*K)] as a Parameter.
    /// May be a ConstantParameter or a MeshElementParameter for
    /// spatially varying values.  Evaluated on-the-fly via SpatialPosition.
    ParameterLib::Parameter<double> const& wall_thermal_conductivity;

    double outsideDiameter() const { return diameter + 2 * wall_thickness; }

    double outsideArea() const { return circleArea(outsideDiameter()); }

    double area() const { return circleArea(diameter); }

    static double circleArea(double const d)
    {
        return std::numbers::pi * d * d / 4;
    }
};

inline double coaxialPipesAnnulusDiameter(Pipe const& inner_pipe,
                                          Pipe const& outer_pipe)
{
    return outer_pipe.diameter - inner_pipe.outsideDiameter();
}

Pipe createPipe(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters);
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
