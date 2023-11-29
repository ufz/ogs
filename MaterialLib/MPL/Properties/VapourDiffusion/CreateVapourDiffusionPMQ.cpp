/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 7, 2021, 9:17 AM
 */

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
