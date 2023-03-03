/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "SaturationVolStrain.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationVolStrain> createSaturationVolStrain(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SaturationVolStrain");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SaturationVolStrain medium property {:s}.", property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__SaturationVolStrain__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__SaturationVolStrain__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const exponent =
        //! \ogs_file_param{properties__property__SaturationVolStrain__exponent}
        config.getConfigParameter<double>("exponent");
    //! \ogs_file_param{properties__property__SaturationVolStrain__p_b}
    auto const p_b = config.getConfigParameter<double>("p_b");
    auto const b11 =
        //! \ogs_file_param{properties__property__SaturationVolStrain__b1}
        config.getConfigParameter<double>("b11");
    auto const b22 =
        //! \ogs_file_param{properties__property__SaturationVolStrain__b2}
        config.getConfigParameter<double>("b22");
    auto const b33 =
        //! \ogs_file_param{properties__property__SaturationVolStrain__b3}
        config.getConfigParameter<double>("b33");

    return std::make_unique<SaturationVolStrain>(
        std::move(property_name), residual_liquid_saturation,
        residual_gas_saturation, exponent, p_b,b11,b22,b33);
}
}  // namespace MaterialPropertyLib
