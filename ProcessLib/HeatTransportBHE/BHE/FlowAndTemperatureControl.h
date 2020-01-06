/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <variant>
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct FlowAndTemperature
{
    double const flow_rate;
    double const temperature;
};

struct TemperatureCurveConstantFlow
{
    FlowAndTemperature operator()(double const /*T_out*/,
                                  double const time) const
    {
        return {flow_rate, temperature_curve.getValue(time)};
    }
    double flow_rate;
    MathLib::PiecewiseLinearInterpolation const& temperature_curve;
};

struct FixedPowerConstantFlow
{
    FlowAndTemperature operator()(double const T_out,
                                  double const /*time*/) const
    {
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    double flow_rate;
    double power;  // Value is expected to be in Watt.
    double heat_capacity;
    double density;
};

struct FixedPowerFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const flow_rate = flow_curve.getValue(time);
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double power;  // Value is expected to be in Watt.
    double heat_capacity;
    double density;
};

struct PowerCurveConstantFlow
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power = power_curve.getValue(time);
        if (power == 0)
        {
            return {0.0, T_out};
        }
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    MathLib::PiecewiseLinearInterpolation const& power_curve;

    double flow_rate;
    double heat_capacity;
    double density;
};

struct BuildingPowerCurveConstantFlow
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power = power_curve.getValue(time);
        if (power == 0)
        {
            return {0.0, T_out};
        }
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    MathLib::PiecewiseLinearInterpolation const& power_curve;

    double flow_rate;
    double heat_capacity;
    double density;
};

using FlowAndTemperatureControl = std::variant<TemperatureCurveConstantFlow,
                                               FixedPowerConstantFlow,
                                               FixedPowerFlowCurve,
                                               PowerCurveConstantFlow,
                                               BuildingPowerCurveConstantFlow>;
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
