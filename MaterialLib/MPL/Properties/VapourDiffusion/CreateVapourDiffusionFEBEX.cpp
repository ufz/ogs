// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateVapourDiffusionFEBEX.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "VapourDiffusionFEBEX.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createVapourDiffusionFEBEX(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "VapourDiffusionFEBEX");
    DBUG("Create VapourDiffusionFEBEX phase property");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    double const base_diffusion_coefficient =
        //! \ogs_file_param{properties__property__VapourDiffusionFEBEX__base_diffusion_coefficient}
        config.getConfigParameter<double>("base_diffusion_coefficient",
                                          2.16e-5);

    double const exponent =
        //! \ogs_file_param{properties__property__VapourDiffusionFEBEX__exponent}
        config.getConfigParameter<double>("exponent", 1.8);

    return std::make_unique<VapourDiffusionFEBEX>(
        std::move(property_name), base_diffusion_coefficient, exponent);
}
}  // namespace MaterialPropertyLib
