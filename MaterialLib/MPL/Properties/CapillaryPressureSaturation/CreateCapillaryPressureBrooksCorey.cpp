/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 25, 2020, 5:00 PM
 */

#include "CreateCapillaryPressureBrooksCorey.h"

#include "BaseLib/ConfigTree.h"
#include "CapillaryPressureBrooksCorey.h"

namespace MaterialPropertyLib
{
std::unique_ptr<CapillaryPressureBrooksCorey>
createCapillaryPressureBrooksCorey(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "CapillaryPressureBrooksCorey");

    DBUG("Create CapillaryPressureBrooksCorey medium property");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__CapillaryPressureBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const maximum_liquid_saturation =
        //! \ogs_file_param{properties__property__CapillaryPressureBrooksCorey__maximum_liquid_saturation}
        config.getConfigParameter<double>("maximum_liquid_saturation");
    auto const exponent =
        //! \ogs_file_param{properties__property__CapillaryPressureBrooksCorey__lambda}
        config.getConfigParameter<double>("lambda");
    auto const entry_pressure =
        //! \ogs_file_param{properties__property__CapillaryPressureBrooksCorey__entry_pressure}
        config.getConfigParameter<double>("entry_pressure");

    auto const maximum_capillary_pressure =
        //! \ogs_file_param{properties__property__CapillaryPressureBrooksCorey__maximum_capillary_pressure}
        config.getConfigParameter<double>("maximum_capillary_pressure");

    return std::make_unique<CapillaryPressureBrooksCorey>(
        residual_liquid_saturation, maximum_liquid_saturation, exponent,
        entry_pressure, maximum_capillary_pressure);
}
}  // namespace MaterialPropertyLib
