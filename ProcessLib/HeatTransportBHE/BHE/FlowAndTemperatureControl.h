// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <variant>

#include "BaseLib/Error.h"
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

inline FlowAndTemperature check_power_and_flow_rate(double flow_rate,
                                                    double power,
                                                    double heat_capacity,
                                                    double density,
                                                    double T_out,
                                                    double flow_rate_min,
                                                    double power_min)
{
    flow_rate = (std::abs(flow_rate) < flow_rate_min) ? 0.0 : flow_rate;

    if (std::abs(power) < power_min)
    {
        return {flow_rate, T_out};
    }

    if (flow_rate == 0)
    {
        OGS_FATAL(
            "BHE power {:0.2f} W with flow rate of 0 m3/s requested, "
            "calculation not possible!",
            power);
    }
    return {flow_rate, power / flow_rate / heat_capacity / density + T_out};
};

struct TemperatureCurveConstantFlow
{
    FlowAndTemperature operator()(double const /*T_out*/,
                                  double const time) const
    {
        return {(std::abs(flow_rate) < flow_rate_min) ? 0.0 : flow_rate,
                temperature_curve.getValue(time)};
    }
    double flow_rate;
    MathLib::PiecewiseLinearInterpolation const& temperature_curve;
    double flow_rate_min;
    static constexpr bool is_power_bc = false;
};

struct TemperatureCurveFlowCurve
{
    FlowAndTemperature operator()(double const /*T_out*/,
                                  double const time) const
    {
        double flow_rate = flow_rate_curve.getValue(time);
        flow_rate = (std::abs(flow_rate) < flow_rate_min) ? 0.0 : flow_rate;
        return {flow_rate, temperature_curve.getValue(time)};
    }
    MathLib::PiecewiseLinearInterpolation const& flow_rate_curve;
    MathLib::PiecewiseLinearInterpolation const& temperature_curve;
    double flow_rate_min;
    static constexpr bool is_power_bc = false;
};

struct FixedPowerConstantFlow
{
    FlowAndTemperature operator()(double const T_out,
                                  double const /*time*/) const
    {
        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    double flow_rate;
    double power;  // Value is expected to be in Watt.
    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct FixedPowerFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const flow_rate = flow_curve.getValue(time);

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double power;  // Value is expected to be in Watt.
    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct PowerCurveConstantFlow
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power = power_curve.getValue(time);

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    MathLib::PiecewiseLinearInterpolation const& power_curve;

    double flow_rate;
    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct PowerCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power = power_curve.getValue(time);
        double const flow_rate = flow_curve.getValue(time);

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    MathLib::PiecewiseLinearInterpolation const& power_curve;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveConstantFlow
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power_building =
            building_power_curves.power_curve.getValue(time);
        double const cop = building_power_curves.cop_curve.getValue(T_out);

        double const power = power_building - power_building / cop;

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    BuildingPowerCurves const building_power_curves;

    double flow_rate;
    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveHotWaterCurveActiveCoolingCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double const power_hot_water =
            building_hot_water_curves.power_curve.getValue(time);
        double const cop_hot_water =
            building_hot_water_curves.cop_curve.getValue(T_out);

        double const power_cooling =
            building_active_cooling_curves.power_curve.getValue(time);
        double const cop_cooling =
            building_active_cooling_curves.cop_curve.getValue(T_out);

        double const flow_rate = flow_curve.getValue(time);

        double const power = power_heating - power_heating / cop_heating +
                             power_hot_water - power_hot_water / cop_hot_water +
                             power_cooling - power_cooling / cop_cooling;

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    BuildingPowerCurves const building_heating_curves;
    BuildingPowerCurves const building_hot_water_curves;
    BuildingPowerCurves const building_active_cooling_curves;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveHotWaterCurvePassiveCoolingCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double const power_hot_water =
            building_hot_water_curves.power_curve.getValue(time);
        double const cop_hot_water =
            building_hot_water_curves.cop_curve.getValue(T_out);

        double const power_cooling = cooling_power.getValue(time);

        double const flow_rate = flow_curve.getValue(time);

        double const power = power_heating - power_heating / cop_heating +
                             power_hot_water - power_hot_water / cop_hot_water +
                             power_cooling;

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    BuildingPowerCurves const building_heating_curves;
    BuildingPowerCurves const building_hot_water_curves;
    MathLib::PiecewiseLinearInterpolation const& cooling_power;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveHotWaterCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double const power_hot_water =
            building_hot_water_curves.power_curve.getValue(time);
        double const cop_hot_water =
            building_hot_water_curves.cop_curve.getValue(T_out);

        double const flow_rate = flow_curve.getValue(time);

        double const power = power_heating - power_heating / cop_heating +
                             power_hot_water - power_hot_water / cop_hot_water;

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    BuildingPowerCurves const building_heating_curves;
    BuildingPowerCurves const building_hot_water_curves;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveActiveCoolingCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double const power_cooling =
            building_active_cooling_curves.power_curve.getValue(time);
        double const cop_cooling =
            building_active_cooling_curves.cop_curve.getValue(T_out);

        double const flow_rate = flow_curve.getValue(time);

        double const power = power_heating - power_heating / cop_heating +
                             power_cooling - power_cooling / cop_cooling;

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    BuildingPowerCurves const building_heating_curves;
    BuildingPowerCurves const building_active_cooling_curves;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurvePassiveCoolingCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double const power_cooling = cooling_power.getValue(time);

        double const flow_rate = flow_curve.getValue(time);

        double const power =
            power_heating - power_heating / cop_heating + power_cooling;

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    BuildingPowerCurves const building_heating_curves;
    MathLib::PiecewiseLinearInterpolation const& cooling_power;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power_heating =
            building_heating_curves.power_curve.getValue(time);
        double const cop_heating =
            building_heating_curves.cop_curve.getValue(T_out);

        double const flow_rate = flow_curve.getValue(time);

        double const power = power_heating - power_heating / cop_heating;

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    BuildingPowerCurves const building_heating_curves;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct ActiveCoolingCurveFlowCurve
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        double const power_cooling =
            building_active_cooling_curves.power_curve.getValue(time);
        double const cop_cooling =
            building_active_cooling_curves.cop_curve.getValue(T_out);

        double const flow_rate = flow_curve.getValue(time);

        double const power = power_cooling - power_cooling / cop_cooling;

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }
    BuildingPowerCurves const building_active_cooling_curves;
    MathLib::PiecewiseLinearInterpolation const& flow_curve;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
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
