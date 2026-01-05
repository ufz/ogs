// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace MathLib
{
class PiecewiseLinearInterpolation;
}

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct BuildingPowerCurves
{
    MathLib::PiecewiseLinearInterpolation const& power_curve;
    MathLib::PiecewiseLinearInterpolation const& cop_curve;
};
using CoolingVariant =
    std::variant<BuildingPowerCurves,
                 std::reference_wrapper<MathLib::PiecewiseLinearInterpolation>>;
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
