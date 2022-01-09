/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "RefrigerantProperties.h"

#include "BaseLib/ConfigTree.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
RefrigerantProperties createRefrigerantProperties(
    BaseLib::ConfigTree const& config)
{
    auto const refrigerant_density =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant__density}
        config.getConfigParameter<double>("density");
    auto const refrigerant_viscosity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant__viscosity}
        config.getConfigParameter<double>("viscosity");
    auto const refrigerant_heat_capacity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant__specific_heat_capacity}
        config.getConfigParameter<double>("specific_heat_capacity");
    auto const refrigerant_thermal_conductivity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant__thermal_conductivity}
        config.getConfigParameter<double>("thermal_conductivity");
    auto const refrigerant_reference_temperature =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__refrigerant__reference_temperature}
        config.getConfigParameter<double>("reference_temperature");
    return {refrigerant_viscosity, refrigerant_density,
            refrigerant_thermal_conductivity, refrigerant_heat_capacity,
            refrigerant_reference_temperature};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
