// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <variant>

#include "BaseLib/Error.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ParameterLib/Parameter.h"
#include "PowerWithCOP.h"

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

struct InflowTemperature
{
    FlowAndTemperature operator()(double const /*T_out*/,
                                  double const time) const
    {
        ParameterLib::SpatialPosition x;
        double flow_rate = flow_rate_param(time, x)[0];
        double temperature = temperature_param(time, x)[0];
        return {(std::abs(flow_rate) < flow_rate_min) ? 0.0 : flow_rate,
                temperature};
    }
    ParameterLib::Parameter<double> const& temperature_param;
    ParameterLib::Parameter<double> const& flow_rate_param;
    double flow_rate_min;
    static constexpr bool is_power_bc = false;
};

struct Power
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        ParameterLib::SpatialPosition x;
        double flow_rate = flow_rate_param(time, x)[0];
        double power = power_param(time, x)[0];
        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }

    // Value is expected to be in Watt.
    ParameterLib::Parameter<double> const& power_param;
    ParameterLib::Parameter<double> const& flow_rate_param;
    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPower
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        ParameterLib::SpatialPosition x;
        double flow_rate = flow_rate_param(time, x)[0];
        double const power_building = building_power.power_param(time, x)[0];
        double const cop = building_power.cop_curve.getValue(T_out);

        double const power = power_building - power_building / cop;

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }

    PowerWithCOP const building_power;
    ParameterLib::Parameter<double> const& flow_rate_param;
    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerHotWaterActiveCooling
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        ParameterLib::SpatialPosition x;
        double flow_rate = flow_rate_param(time, x)[0];
        double const power_heating = building_heating.power_param(time, x)[0];
        double const cop_heating = building_heating.cop_curve.getValue(T_out);

        double const power_hot_water =
            building_hot_water.power_param(time, x)[0];
        double const cop_hot_water =
            building_hot_water.cop_curve.getValue(T_out);

        double const power_cooling =
            building_active_cooling.power_param(time, x)[0];
        double const cop_cooling =
            building_active_cooling.cop_curve.getValue(T_out);

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

    PowerWithCOP const building_heating;
    PowerWithCOP const building_hot_water;
    PowerWithCOP const building_active_cooling;
    ParameterLib::Parameter<double> const& flow_rate_param;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerHotWaterPassiveCooling
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        ParameterLib::SpatialPosition x;
        double flow_rate = flow_rate_param(time, x)[0];
        double const power_heating = building_heating.power_param(time, x)[0];
        double const cop_heating = building_heating.cop_curve.getValue(T_out);

        double const power_hot_water =
            building_hot_water.power_param(time, x)[0];
        double const cop_hot_water =
            building_hot_water.cop_curve.getValue(T_out);

        double const power_cooling = cooling_power_param(time, x)[0];

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

    PowerWithCOP const building_heating;
    PowerWithCOP const building_hot_water;
    ParameterLib::Parameter<double> const& cooling_power_param;
    ParameterLib::Parameter<double> const& flow_rate_param;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerHotWater
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        ParameterLib::SpatialPosition x;
        double flow_rate = flow_rate_param(time, x)[0];
        double const power_heating = building_heating.power_param(time, x)[0];
        double const cop_heating = building_heating.cop_curve.getValue(T_out);

        double const power_hot_water =
            building_hot_water.power_param(time, x)[0];
        double const cop_hot_water =
            building_hot_water.cop_curve.getValue(T_out);

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

    PowerWithCOP const building_heating;
    PowerWithCOP const building_hot_water;
    ParameterLib::Parameter<double> const& flow_rate_param;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerActiveCooling
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        ParameterLib::SpatialPosition x;
        double flow_rate = flow_rate_param(time, x)[0];
        double const power_heating = building_heating.power_param(time, x)[0];
        double const cop_heating = building_heating.cop_curve.getValue(T_out);

        double const power_cooling =
            building_active_cooling.power_param(time, x)[0];
        double const cop_cooling =
            building_active_cooling.cop_curve.getValue(T_out);

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

    PowerWithCOP const building_heating;
    PowerWithCOP const building_active_cooling;
    ParameterLib::Parameter<double> const& flow_rate_param;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct BuildingPowerPassiveCooling
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        ParameterLib::SpatialPosition x;
        double flow_rate = flow_rate_param(time, x)[0];
        double const power_heating = building_heating.power_param(time, x)[0];
        double const cop_heating = building_heating.cop_curve.getValue(T_out);

        double const power_cooling = cooling_power_param(time, x)[0];

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

    PowerWithCOP const building_heating;
    ParameterLib::Parameter<double> const& cooling_power_param;
    ParameterLib::Parameter<double> const& flow_rate_param;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

struct ActiveCooling
{
    FlowAndTemperature operator()(double const T_out, double const time) const
    {
        ParameterLib::SpatialPosition x;
        double const flow_rate = flow_rate_param(time, x)[0];

        double const power_cooling =
            building_active_cooling.power_param(time, x)[0];
        double const cop_cooling =
            building_active_cooling.cop_curve.getValue(T_out);

        double const power = power_cooling - power_cooling / cop_cooling;

        return check_power_and_flow_rate(flow_rate,
                                         power,
                                         heat_capacity,
                                         density,
                                         T_out,
                                         flow_rate_min,
                                         power_min);
    }

    PowerWithCOP const building_active_cooling;
    ParameterLib::Parameter<double> const& flow_rate_param;

    double heat_capacity;
    double density;
    double flow_rate_min;
    double power_min;
    static constexpr bool is_power_bc = true;
};

using FlowAndTemperatureControl =
    std::variant<InflowTemperature,
                 Power,
                 BuildingPower,
                 BuildingPowerHotWaterActiveCooling,
                 BuildingPowerHotWaterPassiveCooling,
                 BuildingPowerHotWater,
                 BuildingPowerActiveCooling,
                 BuildingPowerPassiveCooling,
                 ActiveCooling>;
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
