// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateVapourDiffusionPMQ.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "VapourDiffusionPMQ.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createVapourDiffusionPMQ(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "VapourDiffusionPMQ");
    DBUG("Create VapourDiffusionPMQ phase property");

    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    double const base_diffusion_coefficient =
        //! \ogs_file_param{properties__property__VapourDiffusionPMQ__base_diffusion_coefficient}
        config.getConfigParameter<double>("base_diffusion_coefficient",
                                          2.16e-5);

    double const exponent =
        //! \ogs_file_param{properties__property__VapourDiffusionPMQ__exponent}
        config.getConfigParameter<double>("exponent", 1.8);

    //! \ogs_file_param_special{properties__property__VapourDiffusionPMQ}
    return std::make_unique<VapourDiffusionPMQ>(
        std::move(property_name), base_diffusion_coefficient, exponent);
}
}  // namespace MaterialPropertyLib
