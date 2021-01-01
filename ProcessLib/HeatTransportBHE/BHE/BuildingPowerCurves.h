/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
    MathLib::PiecewiseLinearInterpolation const& cop_heating_curve;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
