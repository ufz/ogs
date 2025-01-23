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

#include <variant>

#include "BuildingPowerCurves.h"
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
    static constexpr bool is_power_bc = false;
};

struct TemperatureCurveFlowCurve
{
    FlowAndTemperature operator()(double const /*T_out*/,
                                  double const time) const
    {
        return {flow_rate_curve.getValue(time),
                temperature_curve.getValue(time)};
    }
    MathLib::PiecewiseLinearInterpolation const& flow_rate_curve;
    MathLib::PiecewiseLinearInterpolation const& temperature_curve;
    static constexpr bool is_power_bc = false;
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
    static constexpr bool is_power_bc = true;
};

struct FixedPowerFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double flow_rate = flow_curve.getValue(time);
        flow_rate = (std::abs(flow_rate) < 1e-12) ? 0.0 : flow_rate;
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double power;  // Value is expected to be in Watt.
    double heat_capacity;
    double density;
    static constexpr bool is_power_bc = true;
};

struct PowerCurveConstantFlow
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double power = power_curve.getValue(time);
        if (std::abs(power) < 1e-12)
        {
            return {flow_rate, T_out};
        }
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    MathLib::PiecewiseLinearInterpolation const& power_curve;

    double flow_rate;
    double heat_capacity;
    double density;
    static constexpr bool is_power_bc = true;
};

struct PowerCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double power = power_curve.getValue(time);
        double flow_rate = flow_curve.getValue(time);
        flow_rate = (std::abs(flow_rate) < 1e-12) ? 0.0 : flow_rate;

        if (std::abs(power) < 1e-12)
        {
            return {flow_rate, T_out};
        }
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    MathLib::PiecewiseLinearInterpolation const& power_curve;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveConstantFlow
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double power = building_power_curves.power_curve.getValue(time);
        double const cop = building_power_curves.cop_curve.getValue(T_out);

        if (std::abs(power) < 1e-12)
        {
            return {flow_rate, T_out};
        }
        return {flow_rate,
                power * (cop - 1) / cop / flow_rate / heat_capacity / density +
                    T_out};
    }
    BuildingPowerCurves const building_power_curves;

    double flow_rate;
    double heat_capacity;
    double density;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveHotWaterCurveActiveCoolingCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double power_hot_water =
            building_hot_water_curves.power_curve.getValue(time);
        double const cop_hot_water =
            building_hot_water_curves.cop_curve.getValue(T_out);

        double power_cooling =
            building_active_cooling_curves.power_curve.getValue(time);
        double const cop_cooling =
            building_active_cooling_curves.cop_curve.getValue(T_out);

        double flow_rate = flow_curve.getValue(time);
        flow_rate = (std::abs(flow_rate) < 1e-12) ? 0.0 : flow_rate;

        double power = power_heating - power_heating / cop_heating +
                       power_hot_water - power_hot_water / cop_hot_water +
                       power_cooling - power_cooling / cop_cooling;

        if (std::abs(power) < 1e-12)
        {
            return {flow_rate, T_out};
        }
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    BuildingPowerCurves const building_heating_curves;
    BuildingPowerCurves const building_hot_water_curves;
    BuildingPowerCurves const building_active_cooling_curves;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveHotWaterCurvePassiveCoolingCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double power_hot_water =
            building_hot_water_curves.power_curve.getValue(time);
        double const cop_hot_water =
            building_hot_water_curves.cop_curve.getValue(T_out);

        double power_cooling = cooling_power.getValue(time);

        double flow_rate = flow_curve.getValue(time);
        flow_rate = (std::abs(flow_rate) < 1e-12) ? 0.0 : flow_rate;

        double power = power_heating - power_heating / cop_heating +
                       power_hot_water - power_hot_water / cop_hot_water +
                       power_cooling;

        if (std::abs(power) < 1e-12)
        {
            return {flow_rate, T_out};
        }
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    BuildingPowerCurves const building_heating_curves;
    BuildingPowerCurves const building_hot_water_curves;
    MathLib::PiecewiseLinearInterpolation const& cooling_power;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveHotWaterCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double power_hot_water =
            building_hot_water_curves.power_curve.getValue(time);
        double const cop_hot_water =
            building_hot_water_curves.cop_curve.getValue(T_out);

        double flow_rate = flow_curve.getValue(time);
        flow_rate = (std::abs(flow_rate) < 1e-12) ? 0.0 : flow_rate;

        double power = power_heating - power_heating / cop_heating +
                       power_hot_water - power_hot_water / cop_hot_water;

        if (std::abs(power) < 1e-12)
        {
            return {flow_rate, T_out};
        }
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    BuildingPowerCurves const building_heating_curves;
    BuildingPowerCurves const building_hot_water_curves;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveActiveCoolingCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double power_cooling =
            building_active_cooling_curves.power_curve.getValue(time);
        double const cop_cooling =
            building_active_cooling_curves.cop_curve.getValue(T_out);

        double flow_rate = flow_curve.getValue(time);
        flow_rate = (std::abs(flow_rate) < 1e-12) ? 0.0 : flow_rate;

        double power = power_heating - power_heating / cop_heating +
                       power_cooling - power_cooling / cop_cooling;

        if (std::abs(power) < 1e-12)
        {
            return {flow_rate, T_out};
        }
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    BuildingPowerCurves const building_heating_curves;
    BuildingPowerCurves const building_active_cooling_curves;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurvePassiveCoolingCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double power_cooling = cooling_power.getValue(time);

        double flow_rate = flow_curve.getValue(time);
        flow_rate = (std::abs(flow_rate) < 1e-12) ? 0.0 : flow_rate;

        double power =
            power_heating - power_heating / cop_heating + power_cooling;

        if (std::abs(power) < 1e-12)
        {
            return {flow_rate, T_out};
        }
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    BuildingPowerCurves const building_heating_curves;
    MathLib::PiecewiseLinearInterpolation const& cooling_power;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double flow_rate = flow_curve.getValue(time);
        flow_rate = (std::abs(flow_rate) < 1e-12) ? 0.0 : flow_rate;

        double power = power_heating - power_heating / cop_heating;

        if (std::abs(power) < 1e-12)
        {
            return {flow_rate, T_out};
        }
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    BuildingPowerCurves const building_heating_curves;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    static constexpr bool is_power_bc = true;
};

struct ActiveCoolingCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double power_cooling =
            building_active_cooling_curves.power_curve.getValue(time);
        double const cop_cooling =
            building_active_cooling_curves.cop_curve.getValue(T_out);

        double flow_rate = flow_curve.getValue(time);
        flow_rate = (std::abs(flow_rate) < 1e-12) ? 0.0 : flow_rate;

        double power = power_cooling - power_cooling / cop_cooling;

        if (std::abs(power) < 1e-12)
        {
            return {flow_rate, T_out};
        }
        return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
    }
    BuildingPowerCurves const building_active_cooling_curves;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    static constexpr bool is_power_bc = true;
};

using FlowAndTemperatureControl =
    std::variant<TemperatureCurveConstantFlow,
                 TemperatureCurveFlowCurve,
                 FixedPowerConstantFlow,
                 FixedPowerFlowCurve,
                 PowerCurveConstantFlow,
                 PowerCurveFlowCurve,
                 BuildingPowerCurveConstantFlow,
                 BuildingPowerCurveHotWaterCurveActiveCoolingCurveFlowCurve,
                 BuildingPowerCurveHotWaterCurvePassiveCoolingCurveFlowCurve,
                 BuildingPowerCurveHotWaterCurveFlowCurve,
                 BuildingPowerCurveActiveCoolingCurveFlowCurve,
                 BuildingPowerCurvePassiveCoolingCurveFlowCurve,
                 BuildingPowerCurveFlowCurve,
                 ActiveCoolingCurveFlowCurve>;
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
