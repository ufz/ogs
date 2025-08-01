/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
