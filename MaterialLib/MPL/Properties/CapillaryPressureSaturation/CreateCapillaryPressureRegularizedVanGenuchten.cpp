/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 21, 2020, 11:06 AM
 */

#include "CreateCapillaryPressureRegularizedVanGenuchten.h"

#include "BaseLib/ConfigTree.h"
#include "CapillaryPressureRegularizedVanGenuchten.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createCapillaryPressureRegularizedVanGenuchten(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type",
                                "CapillaryPressureRegularizedVanGenuchten");

    DBUG("Create CapillaryPressureRegularizedVanGenuchten medium property");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{properties__property__CapillaryPressureRegularizedVanGenuchten__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const maximum_liquid_saturation =
        1.0 -
        //! \ogs_file_param{properties__property__CapillaryPressureRegularizedVanGenuchten__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const exponent =
        //! \ogs_file_param{properties__property__CapillaryPressureRegularizedVanGenuchten__exponent}
        config.getConfigParameter<double>("exponent");
    auto const p_b =
        //! \ogs_file_param{properties__property__CapillaryPressureRegularizedVanGenuchten__p_b}
        config.getConfigParameter<double>("p_b");

    return std::make_unique<CapillaryPressureRegularizedVanGenuchten>(
        residual_liquid_saturation, maximum_liquid_saturation, exponent, p_b);
}
}  // namespace MaterialPropertyLib
