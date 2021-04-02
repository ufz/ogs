/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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

    if (type == "BuildingPowerCurveConstantFlow")
    {
        auto const& power_curve = *BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__BuildingPowerCurveConstantFlow__power_curve}
            config.getConfigParameter<std::string>("power_curve"),
            "Required power curve not found.");

        auto const& cop_heating_curve = *BaseLib::getOrError(
            curves,
            //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__BuildingPowerCurveConstantFlow__cop_heating_curve}
            config.getConfigParameter<std::string>("cop_heating_curve"),
            "Required power curve not found.");

        BuildingPowerCurves const building_power_curves{power_curve,
                                                        cop_heating_curve};

        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__BuildingPowerCurveConstantFlow__flow_rate}
        auto const flow_rate = config.getConfigParameter<double>("flow_rate");

        return BuildingPowerCurveConstantFlow{
            building_power_curves, flow_rate,
            refrigerant.specific_heat_capacity, refrigerant.density};
    }
    OGS_FATAL("FlowAndTemperatureControl type '{:s}' is not implemented.",
              type);
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
