// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
