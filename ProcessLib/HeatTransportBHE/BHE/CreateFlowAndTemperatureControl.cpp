// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateFlowAndTemperatureControl.h"

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "PowerWithCOP.h"
#include "RefrigerantProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
PowerWithCOP createPowerWithCOPStruct(
    std::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPower__hot_water__power}
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPower__heating__power}
    auto const raw_power = config->getConfigParameter<std::string>("power");

    auto const& power_param = ParameterLib::getNamedOrCreateInlineParameter(
        raw_power, parameters, "power", "inline");

    auto const& cop_curve = *BaseLib::getOrError(
        curves,
        //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPower__cooling__cop_curve}
        //! \ogs_file_param_special{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPower__hot_water__cop_curve}
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPower__heating__cop_curve}
        config->getConfigParameter<std::string>("cop_curve"),
        "Required cop curve not found.");

    return PowerWithCOP{power_param, cop_curve};
};

CoolingVariant createCoolingVariant(
    std::optional<BaseLib::ConfigTree> const& cooling_config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    if (  //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPower__cooling__active}
        cooling_config->getConfigParameter<bool>("active", false))
    {
        return createPowerWithCOPStruct(cooling_config, parameters, curves);
    }
    else
    {
        auto const raw_power =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPower__cooling__power}
            cooling_config->getConfigParameter<std::string>("power");

        return ParameterLib::getNamedOrCreateInlineParameter(
            raw_power, parameters, "power", "inline");
    }
};

FlowAndTemperatureControl createHeatingHotWaterCooling(
    std::optional<PowerWithCOP> const& heating,
    std::optional<PowerWithCOP> const& hot_water,
    std::optional<CoolingVariant> const& cooling,
    ParameterLib::Parameter<double> const& flow_rate_param,
    RefrigerantProperties const& refrigerant,
    double const flow_rate_min,
    double const power_min)
{
    if (std::holds_alternative<PowerWithCOP>(*cooling))
    {
        return FlowAndTemperatureControl{
            std::in_place_type<BuildingPowerHotWaterActiveCooling>,
            *heating,
            *hot_water,
            std::get<PowerWithCOP>(*cooling),
            flow_rate_param,
            refrigerant.specific_heat_capacity,
            refrigerant.density,
            flow_rate_min,
            power_min};
    }
    else
    {
        return FlowAndTemperatureControl{
            std::in_place_type<BuildingPowerHotWaterPassiveCooling>,
            *heating,
            *hot_water,
            std::get<std::reference_wrapper<ParameterLib::Parameter<double>>>(
                *cooling)
                .get(),
            flow_rate_param,
            refrigerant.specific_heat_capacity,
            refrigerant.density,
            flow_rate_min,
            power_min};
    }
};

FlowAndTemperatureControl createHeatingCooling(
    std::optional<PowerWithCOP> const& heating,
    std::optional<PowerWithCOP> const& /*hot_water*/,
    std::optional<CoolingVariant> const& cooling,
    ParameterLib::Parameter<double> const& flow_rate_param,
    RefrigerantProperties const& refrigerant,
    double const flow_rate_min,
    double const power_min)
{
    if (std::holds_alternative<PowerWithCOP>(*cooling))
    {
        return FlowAndTemperatureControl{
            std::in_place_type<BuildingPowerActiveCooling>,
            *heating,
            std::get<PowerWithCOP>(*cooling),
            flow_rate_param,
            refrigerant.specific_heat_capacity,
            refrigerant.density,
            flow_rate_min,
            power_min};
    }
    else
    {
        return FlowAndTemperatureControl{
            std::in_place_type<BuildingPowerPassiveCooling>,
            *heating,
            std::get<std::reference_wrapper<ParameterLib::Parameter<double>>>(
                *cooling)
                .get(),
            flow_rate_param,
            refrigerant.specific_heat_capacity,
            refrigerant.density,
            flow_rate_min,
            power_min};
    }
};

FlowAndTemperatureControl createHotWaterCooling(
    std::optional<PowerWithCOP> const& /*heating*/,
    std::optional<PowerWithCOP> const& hot_water,
    std::optional<CoolingVariant> const& cooling,
    ParameterLib::Parameter<double> const& flow_rate_param,
    RefrigerantProperties const& refrigerant,
    double const flow_rate_min,
    double const power_min)
{
    if (std::holds_alternative<PowerWithCOP>(*cooling))
    {
        return FlowAndTemperatureControl{
            std::in_place_type<BuildingPowerActiveCooling>,
            *hot_water,
            std::get<PowerWithCOP>(*cooling),
            flow_rate_param,
            refrigerant.specific_heat_capacity,
            refrigerant.density,
            flow_rate_min,
            power_min};
    }
    else
    {
        return FlowAndTemperatureControl{
            std::in_place_type<BuildingPowerPassiveCooling>,
            *hot_water,
            std::get<std::reference_wrapper<ParameterLib::Parameter<double>>>(
                *cooling)
                .get(),
            flow_rate_param,
            refrigerant.specific_heat_capacity,
            refrigerant.density,
            flow_rate_min,
            power_min};
    }
};

FlowAndTemperatureControl createCooling(
    std::optional<PowerWithCOP> const& /*heating*/,
    std::optional<PowerWithCOP> const& /*hot_water*/,
    std::optional<CoolingVariant> const& cooling,
    ParameterLib::Parameter<double> const& flow_rate_param,
    RefrigerantProperties const& refrigerant,
    double const flow_rate_min,
    double const power_min)
{
    if (std::holds_alternative<PowerWithCOP>(*cooling))
    {
        return FlowAndTemperatureControl{std::in_place_type<ActiveCooling>,
                                         std::get<PowerWithCOP>(*cooling),
                                         flow_rate_param,
                                         refrigerant.specific_heat_capacity,
                                         refrigerant.density,
                                         flow_rate_min,
                                         power_min};
    }
    else
    {
        return FlowAndTemperatureControl{
            std::in_place_type<Power>,
            std::get<std::reference_wrapper<ParameterLib::Parameter<double>>>(
                *cooling)
                .get(),
            flow_rate_param,
            refrigerant.specific_heat_capacity,
            refrigerant.density,
            flow_rate_min,
            power_min};
    }
};

using FactoryAdvancedBuildingCurvesFlowCurve =
    std::function<FlowAndTemperatureControl(
        std::optional<PowerWithCOP>,             // heating
        std::optional<PowerWithCOP>,             // hot water
        std::optional<CoolingVariant>,           // cooling
        ParameterLib::Parameter<double> const&,  // flow rate param
        RefrigerantProperties const&,
        double const,    // flow rate min
        double const)>;  // power min

const std::map<std::tuple<bool, bool, bool>,  // heating, hot_water,
                                              // cooling
               FactoryAdvancedBuildingCurvesFlowCurve>
    advancedBuildingPower = {
        {{true, true, true}, &createHeatingHotWaterCooling},
        {{true, true, false},
         [](std::optional<PowerWithCOP> const& heating,
            std::optional<PowerWithCOP> const& hot_water,
            std::optional<CoolingVariant> const& /*cooling*/,
            auto const& flow_rate_param, auto const& refrigerant,
            auto const flow_rate_min, auto const power_min)
         {
             return FlowAndTemperatureControl{
                 std::in_place_type<BuildingPowerHotWater>,
                 *heating,
                 *hot_water,
                 flow_rate_param,
                 refrigerant.specific_heat_capacity,
                 refrigerant.density,
                 flow_rate_min,
                 power_min};
         }},
        {{true, false, true}, &createHeatingCooling},
        {{false, true, true}, &createHotWaterCooling},
        {{true, false, false},
         [](std::optional<PowerWithCOP> const& heating,
            std::optional<PowerWithCOP> const& /*hot_water*/,
            std::optional<CoolingVariant> const& /*cooling*/,
            auto const& flow_rate_param, auto const& refrigerant,
            auto const flow_rate_min, auto const power_min)
         {
             return FlowAndTemperatureControl{
                 std::in_place_type<BuildingPower>,
                 *heating,
                 flow_rate_param,
                 refrigerant.specific_heat_capacity,
                 refrigerant.density,
                 flow_rate_min,
                 power_min};
         }},
        {{false, true, false},
         [](std::optional<PowerWithCOP> const& /*heating*/,
            std::optional<PowerWithCOP> const& hot_water,
            std::optional<CoolingVariant> const& /*cooling*/,
            auto const& flow_rate_param, auto const& refrigerant,
            auto const flow_rate_min, auto const power_min)
         {
             return FlowAndTemperatureControl{
                 std::in_place_type<BuildingPower>,
                 *hot_water,
                 flow_rate_param,
                 refrigerant.specific_heat_capacity,
                 refrigerant.density,
                 flow_rate_min,
                 power_min};
         }},
        {{false, false, true}, &createCooling}};

FlowAndTemperatureControl createAdvancedBuildingPower(
    std::optional<BaseLib::ConfigTree> const& heating_config,
    std::optional<BaseLib::ConfigTree> const& hot_water_config,
    std::optional<BaseLib::ConfigTree> const& cooling_config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    ParameterLib::Parameter<double> const& flow_rate_param,
    RefrigerantProperties const& refrigerant,
    double const flow_rate_min,
    double const power_min)
{
    std::optional<PowerWithCOP> building_heating_curves;
    std::optional<PowerWithCOP> building_hot_water_curves;
    std::optional<CoolingVariant> building_cooling_curves;

    bool heating = false;
    bool hot_water = false;
    bool cooling = false;

    if (heating_config)
    {
        building_heating_curves.emplace(
            createPowerWithCOPStruct(heating_config, parameters, curves));
        heating = true;
    }
    if (hot_water_config)
    {
        building_hot_water_curves.emplace(
            createPowerWithCOPStruct(hot_water_config, parameters, curves));
        hot_water = true;
    }
    if (cooling_config)
    {
        building_cooling_curves.emplace(
            createCoolingVariant(cooling_config, parameters, curves));
        cooling = true;
    }
    auto key = std::make_tuple(heating, hot_water, cooling);

    auto it = advancedBuildingPower.find(key);
    if (it == advancedBuildingPower.end())
        OGS_FATAL(
            "AdvancedBuildingPower combination is not "
            "implemented.");
    auto factory = it->second;

    return factory(building_heating_curves,
                   building_hot_water_curves,
                   building_cooling_curves,
                   flow_rate_param,
                   refrigerant,
                   flow_rate_min,
                   power_min);
}
FlowAndTemperatureControl createFlowAndTemperatureControl(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    RefrigerantProperties const& refrigerant)
{
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__type}
    auto const type = config.getConfigParameter<std::string>("type");

    auto const flow_rate_min =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__flow_rate_min}
        config.getConfigParameter<double>("flow_rate_min", 1e-6);
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__power_min}
    auto const power_min = config.getConfigParameter<double>("power_min", 1e-3);
    if (type == "InflowTemperature")
    {
        auto const raw_temperature =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__InflowTemperature__temperature}
            config.getConfigParameter<std::string>("temperature");

        auto const& temperature_param =
            ParameterLib::getNamedOrCreateInlineParameter(
                raw_temperature, parameters, "temperature", "inline");

        auto const raw_flow_rate =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__InflowTemperature__flow_rate}
            config.getConfigParameter<std::string>("flow_rate");

        auto const& flow_rate_param =
            ParameterLib::getNamedOrCreateInlineParameter(
                raw_flow_rate, parameters, "flow_rate", "inline");

        return InflowTemperature{temperature_param, flow_rate_param,
                                 flow_rate_min};
    }
    if (type == "Power")
    {
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__Power__power}
        auto const raw_power = config.getConfigParameter<std::string>("power");

        auto const& power_param = ParameterLib::getNamedOrCreateInlineParameter(
            raw_power, parameters, "power", "inline");

        auto const raw_flow_rate =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__Power__flow_rate}
            config.getConfigParameter<std::string>("flow_rate");

        auto const& flow_rate_param =
            ParameterLib::getNamedOrCreateInlineParameter(
                raw_flow_rate, parameters, "flow_rate", "inline");

        return Power{power_param,
                     flow_rate_param,
                     refrigerant.specific_heat_capacity,
                     refrigerant.density,
                     flow_rate_min,
                     power_min};
    }

    if (type == "BuildingPower")
    {
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__BuildingPower__power}
        auto const raw_power = config.getConfigParameter<std::string>("power");

        auto const& power_param = ParameterLib::getNamedOrCreateInlineParameter(
            raw_power, parameters, "power", "inline");

        auto const& cop_curve = *BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__BuildingPower__cop_curve}
            config.getConfigParameter<std::string>("cop_curve"),
            "Required cop curve not found.");

        PowerWithCOP const building_power{power_param, cop_curve};

        auto const raw_flow_rate =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__BuildingPower__flow_rate}
            config.getConfigParameter<std::string>("flow_rate");

        auto const& flow_rate_param =
            ParameterLib::getNamedOrCreateInlineParameter(
                raw_flow_rate, parameters, "flow_rate", "inline");

        return BuildingPower{building_power,
                             flow_rate_param,
                             refrigerant.specific_heat_capacity,
                             refrigerant.density,
                             flow_rate_min,
                             power_min};
    }

    if (type == "AdvancedBuildingPower")
    {
        // add a heating, hot water and cooling config as optional to handle
        auto const& heating_config =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPower__heating}
            config.getConfigSubtreeOptional(
                "heating");  // Optional, take care if it is not present

        auto const& hot_water_config =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPower__hot_water}
            config.getConfigSubtreeOptional(
                "hot_water");  // Optional, take care if it is not present

        auto const& cooling_config =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPower__cooling}
            config.getConfigSubtreeOptional(
                "cooling");  // Optional, take care if it is not present

        auto const raw_flow_rate =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPower__flow_rate}
            config.getConfigParameter<std::string>("flow_rate");

        auto const& flow_rate_param =
            ParameterLib::getNamedOrCreateInlineParameter(
                raw_flow_rate, parameters, "flow_rate", "inline");

        return createAdvancedBuildingPower(
            heating_config, hot_water_config, cooling_config, parameters,
            curves, flow_rate_param, refrigerant, flow_rate_min, power_min);
    }
    OGS_FATAL("FlowAndTemperatureControl type '{:s}' is not implemented.",
              type);
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
