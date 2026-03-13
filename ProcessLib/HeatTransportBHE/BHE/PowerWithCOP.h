// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <variant>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct PowerWithCOP
{
    ParameterLib::Parameter<double> const& power_param;
    MathLib::PiecewiseLinearInterpolation const& cop_curve;
};
using CoolingVariant =
    std::variant<PowerWithCOP,
                 std::reference_wrapper<ParameterLib::Parameter<double>>>;
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
