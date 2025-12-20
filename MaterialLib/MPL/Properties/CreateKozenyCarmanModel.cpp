// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateKozenyCarmanModel.h"

#include "BaseLib/ConfigTree.h"
#include "KozenyCarmanModel.h"
#include "ParameterLib/Utils.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createKozenyCarmanModel(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "KozenyCarman");
    DBUG("Create Kozeny-Carman model.");

    auto const& k0 = ParameterLib::findParameter<double>(
        //! \ogs_file_param{properties__property__KozenyCarman__initial_permeability}
        config.getConfigParameter<std::string>("initial_permeability"),
        parameters, 0, nullptr);

    auto const& phi0 = ParameterLib::findParameter<double>(
        //! \ogs_file_param{properties__property__KozenyCarman__initial_porosity}
        config.getConfigParameter<std::string>("initial_porosity"), parameters,
        1, nullptr);

    return std::make_unique<KozenyCarmanModel>(k0, phi0);
}
}  // namespace MaterialPropertyLib
