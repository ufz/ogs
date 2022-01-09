/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BoreholeGeometry.h"

#include "BaseLib/ConfigTree.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
BoreholeGeometry createBoreholeGeometry(BaseLib::ConfigTree const& config)
{
    const auto borehole_length =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole__length}
        config.getConfigParameter<double>("length");
    const auto borehole_diameter =
        //! \ogs_file_param{prj__processes__process__HEAT_TRANSPORT_BHE__borehole_heat_exchangers__borehole_heat_exchanger__borehole__diameter}
        config.getConfigParameter<double>("diameter");
    return {borehole_length, borehole_diameter};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
