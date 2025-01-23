/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateFlowAndTemperatureControl.h"

#include "BaseLib/Algorithm.h"
#include "BaseLib/ConfigTree.h"
#include "BuildingPowerCurves.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "RefrigerantProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
BuildingPowerCurves createBuildingPowerCurvesStruct(
    std::optional<BaseLib::ConfigTree> const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    auto const& power_curve = *BaseLib::getOrError(
        curves,
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPowerCurvesFlowCurve__heating__power_curve}
        config->getConfigParameter<std::string>("power_curve"),
        "Required power curve not found.");

    auto const& cop_curve = *BaseLib::getOrError(
        curves,
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__BuildingPowerCurveConstantFlow__heating__cop_curve}
        config->getConfigParameter<std::string>("cop_curve"),
        "Required cop curve not found.");

    return BuildingPowerCurves{power_curve, cop_curve};
};

CoolingVariant createCoolingVariant(
    std::optional<BaseLib::ConfigTree> const& cooling_config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    if (  //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPowerCurvesFlowCurve__cooling__active}
        cooling_config->getConfigParameter<bool>("active", false))
    {
        return createBuildingPowerCurvesStruct(cooling_config, curves);
    }
    else
    {
        return std::ref(*BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPowerCurvesFlowCurve__cooling__power_curve}
            cooling_config->getConfigParameter<std::string>("power_curve"),
            "Required cooling power curve not found."));
    }
};

FlowAndTemperatureControl createHeatingHotWaterCooling(
    BuildingPowerCurves const& heating,
    BuildingPowerCurves const& hot_water,
    CoolingVariant const& cooling,
    MathLib::PiecewiseLinearInterpolation const& flow_rate_curve,
    RefrigerantProperties const& refrigerant)
{
    if (std::holds_alternative<BuildingPowerCurves>(cooling))
    {
        return FlowAndTemperatureControl{
            std::in_place_type<
                BuildingPowerCurveHotWaterCurveActiveCoolingCurveFlowCurve>,
            heating,
            hot_water,
            std::get<BuildingPowerCurves>(cooling),
            flow_rate_curve,
            refrigerant.specific_heat_capacity,
            refrigerant.density};
    }
    else
    {
        return FlowAndTemperatureControl{
            std::in_place_type<
                BuildingPowerCurveHotWaterCurvePassiveCoolingCurveFlowCurve>,
            heating,
            hot_water,
            std::get<
                std::reference_wrapper<MathLib::PiecewiseLinearInterpolation>>(
                cooling)
                .get(),
            flow_rate_curve,
            refrigerant.specific_heat_capacity,
            refrigerant.density};
    }
};

FlowAndTemperatureControl createHeatingCooling(
    BuildingPowerCurves const& heating,
    BuildingPowerCurves const& /*hot_water*/,
    CoolingVariant const& cooling,
    MathLib::PiecewiseLinearInterpolation const& flow_rate_curve,
    RefrigerantProperties const& refrigerant)
{
    if (std::holds_alternative<BuildingPowerCurves>(cooling))
    {
        return FlowAndTemperatureControl{
            std::in_place_type<BuildingPowerCurveActiveCoolingCurveFlowCurve>,
            heating,
            std::get<BuildingPowerCurves>(cooling),
            flow_rate_curve,
            refrigerant.specific_heat_capacity,
            refrigerant.density};
    }
    else
    {
        return FlowAndTemperatureControl{
            std::in_place_type<BuildingPowerCurvePassiveCoolingCurveFlowCurve>,
            heating,
            std::get<
                std::reference_wrapper<MathLib::PiecewiseLinearInterpolation>>(
                cooling)
                .get(),
            flow_rate_curve,
            refrigerant.specific_heat_capacity,
            refrigerant.density};
    }
};

FlowAndTemperatureControl createHotWaterCooling(
    BuildingPowerCurves const& /*heating*/,
    BuildingPowerCurves const& hot_water,
    CoolingVariant const& cooling,
    MathLib::PiecewiseLinearInterpolation const& flow_rate_curve,
    RefrigerantProperties const& refrigerant)
{
    if (std::holds_alternative<BuildingPowerCurves>(cooling))
    {
        return FlowAndTemperatureControl{
            std::in_place_type<BuildingPowerCurveActiveCoolingCurveFlowCurve>,
            hot_water,
            std::get<BuildingPowerCurves>(cooling),
            flow_rate_curve,
            refrigerant.specific_heat_capacity,
            refrigerant.density};
    }
    else
    {
        return FlowAndTemperatureControl{
            std::in_place_type<BuildingPowerCurvePassiveCoolingCurveFlowCurve>,
            hot_water,
            std::get<
                std::reference_wrapper<MathLib::PiecewiseLinearInterpolation>>(
                cooling)
                .get(),
            flow_rate_curve,
            refrigerant.specific_heat_capacity,
            refrigerant.density};
    }
};

FlowAndTemperatureControl createCooling(
    BuildingPowerCurves const& /*heating*/,
    BuildingPowerCurves const& /*hot_water*/,
    CoolingVariant const& cooling,
    MathLib::PiecewiseLinearInterpolation const& flow_rate_curve,
    RefrigerantProperties const& refrigerant)
{
    if (std::holds_alternative<BuildingPowerCurves>(cooling))
    {
        return FlowAndTemperatureControl{
            std::in_place_type<ActiveCoolingCurveFlowCurve>,
            std::get<BuildingPowerCurves>(cooling), flow_rate_curve,
            refrigerant.specific_heat_capacity, refrigerant.density};
    }
    else
    {
        return FlowAndTemperatureControl{
            std::in_place_type<PowerCurveFlowCurve>,
            std::get<
                std::reference_wrapper<MathLib::PiecewiseLinearInterpolation>>(
                cooling)
                .get(),
            flow_rate_curve, refrigerant.specific_heat_capacity,
            refrigerant.density};
    }
};

using FactoryAdvancedBuildingCurvesFlowCurve =
    std::function<FlowAndTemperatureControl(
        BuildingPowerCurves const&,                    // heating
        BuildingPowerCurves const&,                    // hot water
        CoolingVariant const&,                         // cooling
        MathLib::PiecewiseLinearInterpolation const&,  // flow rate curve
        RefrigerantProperties const&)>;

const std::map<std::tuple<bool, bool, bool>,  // heating, hot_water,
                                              // cooling
               FactoryAdvancedBuildingCurvesFlowCurve>
    advancedBuildingPowerCurvesFlowCurve = {
        {{true, true, true}, &createHeatingHotWaterCooling},
        {{true, true, false},
         [](auto const& heating, auto const& hot_water, auto const& /*cooling*/,
            auto const& flow_rate_curve, auto const& refrigerant)
         {
             return FlowAndTemperatureControl{
                 std::in_place_type<BuildingPowerCurveHotWaterCurveFlowCurve>,
                 heating,
                 hot_water,
                 flow_rate_curve,
                 refrigerant.specific_heat_capacity,
                 refrigerant.density};
         }},
        {{true, false, true}, &createHeatingCooling},
        {{false, true, true}, &createHotWaterCooling},
        {{true, false, false},
         [](auto const& heating, auto const& /*hot_water*/,
            auto const& /*cooling*/, auto const& flow_rate_curve,
            auto const& refrigerant)
         {
             return FlowAndTemperatureControl{
                 std::in_place_type<BuildingPowerCurveFlowCurve>, heating,
                 flow_rate_curve, refrigerant.specific_heat_capacity,
                 refrigerant.density};
         }},
        {{false, true, false},
         [](auto const& /*heating*/, auto const& hot_water,
            auto const& /*cooling*/, auto const& flow_rate_curve,
            auto const& refrigerant)
         {
             return FlowAndTemperatureControl{
                 std::in_place_type<BuildingPowerCurveFlowCurve>, hot_water,
                 flow_rate_curve, refrigerant.specific_heat_capacity,
                 refrigerant.density};
         }},
        {{false, false, true}, &createCooling}};

FlowAndTemperatureControl createAdvancedBuildingPowerCurvesFlowCurve(
    std::optional<BaseLib::ConfigTree> const& heating_config,
    std::optional<BaseLib::ConfigTree> const& hot_water_config,
    std::optional<BaseLib::ConfigTree> const& cooling_config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    MathLib::PiecewiseLinearInterpolation const& flow_rate_curve,
    RefrigerantProperties const& refrigerant)
{
    std::optional<BuildingPowerCurves> building_heating_curves;
    std::optional<BuildingPowerCurves> building_hot_water_curves;
    std::optional<CoolingVariant> building_cooling_curves;

    bool heating = false;
    bool hot_water = false;
    bool cooling = false;

    if (heating_config)
    {
        building_heating_curves.emplace(
            createBuildingPowerCurvesStruct(heating_config, curves));
        heating = true;
    }
    if (hot_water_config)
    {
        building_hot_water_curves.emplace(
            createBuildingPowerCurvesStruct(hot_water_config, curves));
        hot_water = true;
    }
    if (cooling_config)
    {
        building_cooling_curves.emplace(
            createCoolingVariant(cooling_config, curves));
        cooling = true;
    }
    auto key = std::make_tuple(heating, hot_water, cooling);

    auto it = advancedBuildingPowerCurvesFlowCurve.find(key);
    if (it == advancedBuildingPowerCurvesFlowCurve.end())
        OGS_FATAL(
            "AdvancedBuildingPowerCurvesFlowCurve combination is not "
            "implemented.");
    auto factory = it->second;

    return factory(*building_heating_curves,
                   *building_hot_water_curves,
                   *building_cooling_curves,
                   flow_rate_curve,
                   refrigerant);
}
FlowAndTemperatureControl createFlowAndTemperatureControl(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    RefrigerantProperties const& refrigerant)
{
    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__type}
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "TemperatureCurveConstantFlow")
    {
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__TemperatureCurveConstantFlow__flow_rate}
        auto const flow_rate = config.getConfigParameter<double>("flow_rate");

        auto const& temperature_curve = *BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__TemperatureCurveConstantFlow__temperature_curve}
            config.getConfigParameter<std::string>("temperature_curve"),
            "Required temperature curve not found.");

        return TemperatureCurveConstantFlow{flow_rate, temperature_curve};
    }
    if (type == "TemperatureCurveFlowCurve")
    {
        auto const& flow_rate_curve = *BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__TemperatureCurveFlowCurve__flow_rate_curve}
            config.getConfigParameter<std::string>("flow_rate_curve"),
            "Required flow curve not found.");

        auto const& temperature_curve = *BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__TemperatureCurveFlowCurve__temperature_curve}
            config.getConfigParameter<std::string>("temperature_curve"),
            "Required temperature curve not found.");

        return TemperatureCurveFlowCurve{flow_rate_curve, temperature_curve};
    }
    if (type == "FixedPowerConstantFlow")
    {
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__FixedPowerConstantFlow__power}
        auto const power = config.getConfigParameter<double>("power");

        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__FixedPowerConstantFlow__flow_rate}
        auto const flow_rate = config.getConfigParameter<double>("flow_rate");

        return FixedPowerConstantFlow{flow_rate, power,
                                      refrigerant.specific_heat_capacity,
                                      refrigerant.density};
    }

    if (type == "FixedPowerFlowCurve")
    {
        auto const& flow_rate_curve = *BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__FixedPowerFlowCurve__flow_rate_curve}
            config.getConfigParameter<std::string>("flow_rate_curve"),
            "Required flow rate curve not found.");

        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__FixedPowerFlowCurve__power}
        auto const power = config.getConfigParameter<double>("power");

        return FixedPowerFlowCurve{flow_rate_curve, power,
                                   refrigerant.specific_heat_capacity,
                                   refrigerant.density};
    }

    if (type == "PowerCurveConstantFlow")
    {
        auto const& power_curve = *BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__PowerCurveConstantFlow__power_curve}
            config.getConfigParameter<std::string>("power_curve"),
            "Required power curve not found.");

        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__PowerCurveConstantFlow__flow_rate}
        auto const flow_rate = config.getConfigParameter<double>("flow_rate");

        return PowerCurveConstantFlow{power_curve, flow_rate,
                                      refrigerant.specific_heat_capacity,
                                      refrigerant.density};
    }

    if (type == "PowerCurveFlowCurve")
    {
        auto const& power_curve = *BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__PowerCurveFlowCurve__power_curve}
            config.getConfigParameter<std::string>("power_curve"),
            "Required power curve not found.");

        auto const& flow_rate_curve = *BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__PowerCurveFlowCurve__flow_rate_curve}
            config.getConfigParameter<std::string>("flow_rate_curve"),
            "Required flow rate curve not found.");

        return PowerCurveFlowCurve{power_curve, flow_rate_curve,
                                   refrigerant.specific_heat_capacity,
                                   refrigerant.density};
    }

    if (type == "AdvancedBuildingPowerCurvesFlowCurve")
    {
        // add a heating, hot water and cooling config as optional to handle
        // different combinations later
        auto const& heating_config =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__flow_and_temperature_control__AdvancedBuildingPowerCurvesFlowCurve}
            config.getConfigSubtreeOptional(
                "heating");  // Optional, take care if it is not present

        // add a heating config to differ between different types
        auto const& hot_water_config =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__flow_and_temperature_control__AdvancedBuildingPowerCurvesFlowCurve}
            config.getConfigSubtreeOptional(
                "hot_water");  // Optional, take care if it is not present

        // add a heating config to differ between different types
        auto const& cooling_config =
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__flow_and_temperature_control__AdvancedBuildingPowerCurvesFlowCurve}
            config.getConfigSubtreeOptional(
                "cooling");  // Optional, take care if it is not present

        auto const& flow_rate_curve = *BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__AdvancedBuildingPowerCurvesFlowCurve__flow_rate_curve}
            config.getConfigParameter<std::string>("flow_rate_curve"),
            "Required flow rate curve not found.");

        return createAdvancedBuildingPowerCurvesFlowCurve(heating_config,
                                                          hot_water_config,
                                                          cooling_config,
                                                          curves,
                                                          flow_rate_curve,
                                                          refrigerant);
    }
    OGS_FATAL("FlowAndTemperatureControl type '{:s}' is not implemented.",
              type);
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
