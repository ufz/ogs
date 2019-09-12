/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "RelPermBrooksCorey.h"

namespace MaterialPropertyLib
{
std::unique_ptr<RelPermBrooksCorey> createRelPermBrooksCorey(
    BaseLib::ConfigTree const& config)
{
    // check is reading the parameter, not peeking it...
    //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey}
    config.checkConfigParameter("type", "RelPermBrooksCorey");
    DBUG("Create RelPermBrooksCorey medium property");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const min_relative_permeability_liquid =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__k_rel_min_liquid}
        config.getConfigParameter<double>("min_relative_permeability_liquid");
    auto const min_relative_permeability_gas =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__k_rel_min_gas}
        config.getConfigParameter<double>("min_relative_permeability_gas");
    auto const exponent =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__lambda}
        config.getConfigParameter<double>("lambda");
    if (exponent <= 0.)
    {
        OGS_FATAL("Exponent 'lambda' must be positive.");
    }

    return std::make_unique<RelPermBrooksCorey>(
        residual_liquid_saturation,
        residual_gas_saturation,
        min_relative_permeability_liquid,
        min_relative_permeability_gas,
        exponent);
}
}  // namespace MaterialPropertyLib
