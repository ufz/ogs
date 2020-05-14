/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 20, 2020, 1:31 PM
 */

#include "CreateCapillaryPressureVanGenuchten.h"

#include <limits>

#include "BaseLib/ConfigTree.h"
#include "CapillaryPressureVanGenuchten.h"
#include "GetSaturationVanGenuchten.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createCapillaryPressureVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "CapillaryPressureVanGenuchten");

    DBUG("Create CapillaryPressureVanGenuchten medium property");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const exponent =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__exponent}
        config.getConfigParameter<double>("exponent");
    auto const p_b =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__p_b}
        config.getConfigParameter<double>("p_b");
    auto const maximum_capillary_pressure =
        //! \ogs_file_param{properties__property__CapillaryPressureVanGenuchten__maximum_capillary_pressure}
        config.getConfigParameter<double>(
            "maximum_capillary_pressure",
            std::numeric_limits<double>::quiet_NaN());

    // If maximum_capillary_pressure, pc_max, is not given,
    //    S_at_pc_max = S_L_res  + 1.0e-9.
    // Otherwise
    //    S_at_pc_max = max(S_L_res  + 1.0e-9, S(pc_max))
    const double S_at_pc_max =
        (maximum_capillary_pressure == std::numeric_limits<double>::quiet_NaN())
            ? residual_liquid_saturation + 1.0e-9
            : std::max(residual_liquid_saturation + 1.0e-9,
                       getSaturationVanGenuchten(
                           maximum_capillary_pressure, p_b,
                           residual_liquid_saturation,
                           1.0 - residual_gas_saturation, exponent));

    return std::make_unique<CapillaryPressureVanGenuchten>(
        residual_liquid_saturation, residual_gas_saturation, exponent, p_b,
        S_at_pc_max);
}
}  // namespace MaterialPropertyLib
