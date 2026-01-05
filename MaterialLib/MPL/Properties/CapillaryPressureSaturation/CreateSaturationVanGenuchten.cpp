// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BaseLib/ConfigTree.h"
#include "SaturationVanGenuchten.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationVanGenuchten> createSaturationVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SaturationVanGenuchten");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SaturationVanGenuchten medium property {:s}.", property_name);

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__SaturationVanGenuchten__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__SaturationVanGenuchten__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");

    double pressure_exponent;
    double saturation_exponent;
    if (auto const optional_exponent =
            //! \ogs_file_param{properties__property__SaturationVanGenuchten__exponent}
        config.getConfigParameterOptional<double>("exponent"))
    {
        pressure_exponent = *optional_exponent;
        saturation_exponent = 1.0 / (1.0 - pressure_exponent);
    }
    else
    {
        pressure_exponent =
            //! \ogs_file_param{properties__property__SaturationVanGenuchten__pressure_exponent}
            config.getConfigParameter<double>("pressure_exponent");
        saturation_exponent =
            //! \ogs_file_param{properties__property__SaturationVanGenuchten__saturation_exponent}
            config.getConfigParameter<double>("saturation_exponent");
    }

    //! \ogs_file_param{properties__property__SaturationVanGenuchten__p_b}
    auto const p_b = config.getConfigParameter<double>("p_b");

    return std::make_unique<SaturationVanGenuchten>(
        std::move(property_name), residual_liquid_saturation,
        residual_gas_saturation, pressure_exponent, saturation_exponent, p_b);
}
}  // namespace MaterialPropertyLib
