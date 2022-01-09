/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "RelPermBrooksCorey");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create RelPermBrooksCorey medium property {:s}.", property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__RelPermBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__RelPermBrooksCorey__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const min_relative_permeability =
        //! \ogs_file_param{properties__property__RelPermBrooksCorey__min_relative_permeability}
        config.getConfigParameter<double>("min_relative_permeability");
    auto const exponent =
        //! \ogs_file_param{properties__property__RelPermBrooksCorey__lambda}
        config.getConfigParameter<double>("lambda");
    if (exponent <= 0.)
    {
        OGS_FATAL("Exponent 'lambda' must be positive.");
    }

    return std::make_unique<RelPermBrooksCorey>(std::move(property_name),
                                                residual_liquid_saturation,
                                                residual_gas_saturation,
                                                min_relative_permeability,
                                                exponent);
}
}  // namespace MaterialPropertyLib
