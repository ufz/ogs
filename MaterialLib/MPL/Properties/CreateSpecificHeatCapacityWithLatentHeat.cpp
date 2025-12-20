// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateSpecificHeatCapacityWithLatentHeat.h"

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/Utils.h"
#include "SpecificHeatCapacityWithLatentHeat.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SpecificHeatCapacityWithLatentHeat>
createSpecificHeatCapacityWithLatentHeat(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SpecificHeatCapacityWithLatentHeat");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create temperature dependent specific heat capacity {:s}.",
         property_name);

    auto const l =
        //! \ogs_file_param{properties__property__SpecificHeatCapacityWithLatentHeat__specific_latent_heat}
        config.getConfigParameter<double>("specific_latent_heat");

    return std::make_unique<
        MaterialPropertyLib::SpecificHeatCapacityWithLatentHeat>(
        std::move(property_name), l);
}
}  // namespace MaterialPropertyLib
