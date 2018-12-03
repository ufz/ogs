/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "BaseLib/Algorithm.h"

#include "CreateFlowAndTemperatureControl.h"
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
        double const flow_rate = config.getConfigParameter<double>("flow_rate");

        auto const& temperature_curve = *BaseLib::getOrError(
                    curves,
                    //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__TemperatureCurveConstantFlow__temperature_curve}
                    config.getConfigParameter<std::string>("temperature_curve"),
                    "Required temperature curve not found.");

        return TemperatureCurveConstantFlow{flow_rate, temperature_curve};
    }
    if (type == "FixedPowerConstantFlow")
    {
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__FixedPowerConstantFlow__power}
        double const power = config.getConfigParameter<double>("power");

        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__flow_and_temperature_control__FixedPowerConstantFlow__flow_rate}
        double const flow_rate = config.getConfigParameter<double>("flow_rate");
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
        double const power = config.getConfigParameter<double>("power");

        return FixedPowerFlowCurve{flow_rate_curve, power,
                                   refrigerant.specific_heat_capacity,
                                   refrigerant.density};
    }
    OGS_FATAL("FlowAndTemperatureControl type '%s' is not implemented.",
              type.c_str());
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
