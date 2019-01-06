/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/variant.hpp>
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

using FlowAndTemperatureControl = boost::variant<TemperatureCurveConstantFlow,
                                                 FixedPowerConstantFlow,
                                                 FixedPowerFlowCurve>;

}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
