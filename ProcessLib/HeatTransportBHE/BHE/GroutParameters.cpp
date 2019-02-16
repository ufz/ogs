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

#include "GroutParameters.h"
#include "BaseLib/ConfigTree.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
GroutParameters createGroutParameters(BaseLib::ConfigTree const& config)
{
    const double grout_density =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout__density}
        config.getConfigParameter<double>("density");
    const double grout_porosity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout__porosity}
        config.getConfigParameter<double>("porosity");
    const double grout_heat_capacity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout__heat_capacity}
        config.getConfigParameter<double>("heat_capacity");
    const double grout_thermal_conductivity =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__grout__thermal_conductivity}
        config.getConfigParameter<double>("thermal_conductivity");
    return {grout_density, grout_porosity, grout_heat_capacity,
            grout_thermal_conductivity};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
