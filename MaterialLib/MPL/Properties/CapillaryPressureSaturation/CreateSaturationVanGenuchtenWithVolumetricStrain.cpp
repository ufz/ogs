/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "SaturationVanGenuchtenWithVolumetricStrain.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationVanGenuchtenWithVolumetricStrain>
createSaturationVanGenuchtenWithVolumetricStrain(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type",
                                "SaturationVanGenuchtenWithVolumetricStrain");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG(
        "Create SaturationVanGenuchtenWithVolumetricStrain medium property "
        "{:s}.",
        property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__SaturationVanGenuchtenWithVolumetricStrain__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__SaturationVanGenuchtenWithVolumetricStrain__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const exponent =
        //! \ogs_file_param{properties__property__SaturationVanGenuchtenWithVolumetricStrain__exponent}
        config.getConfigParameter<double>("exponent");
    //! \ogs_file_param{properties__property__SaturationVanGenuchtenWithVolumetricStrain__p_b}
    auto const p_b = config.getConfigParameter<double>("p_b");
    auto const e_0 =
        //! \ogs_file_param{properties__property__SaturationVanGenuchtenWithVolumetricStrain__e_0}
        config.getConfigParameter<double>("e_0");
    auto const e_m =
        //! \ogs_file_param{properties__property__SaturationVanGenuchtenWithVolumetricStrain__e_m}
        config.getConfigParameter<double>("e_m");
    auto const a =
        //! \ogs_file_param{properties__property__SaturationVanGenuchtenWithVolumetricStrain__a}
        config.getConfigParameter<double>("a");
    auto const d_diff =
        //! \ogs_file_param{properties__property__SaturationVanGenuchtenWithVolumetricStrain__d_diff}
        config.getConfigParameter<double>("d_diff");

    return std::make_unique<SaturationVanGenuchtenWithVolumetricStrain>(
        std::move(property_name), residual_liquid_saturation,
        residual_gas_saturation, exponent, p_b, e_0, e_m, a, d_diff);
}
}  // namespace MaterialPropertyLib
