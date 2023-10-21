/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 30, 2023, 1:51 PM
 */

#include "CreateVapourDiffusionDeVries.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "VapourDiffusionDeVries.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createVapourDiffusionDeVries(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "VapourDiffusionDeVries");
    DBUG("Create VapourDiffusionDeVries phase property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    double const base_diffusion_coefficient =
        //! \ogs_file_param{properties__property__VapourDiffusionDeVries__base_diffusion_coefficient}
        config.getConfigParameter<double>("base_diffusion_coefficient", 5.9e-6);

    double const exponent =
        //! \ogs_file_param{properties__property__VapourDiffusionDeVries__exponent}
        config.getConfigParameter<double>("exponent", 2.3);

    return std::make_unique<VapourDiffusionDeVries>(
        std::move(property_name), base_diffusion_coefficient, exponent);
}
}  // namespace MaterialPropertyLib
