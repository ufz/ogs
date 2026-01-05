// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateVermaPruessModel.h"

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/Utils.h"
#include "VermaPruessModel.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createVermaPruessModel(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "VermaPruess");
    DBUG("Create Verma-Pruess model.");

    auto const& k0 = ParameterLib::findParameter<double>(
        //! \ogs_file_param{properties__property__VermaPruessModel__initial_permeability}
        config.getConfigParameter<std::string>("initial_permeability"),
        parameters, 0, nullptr);

    auto const& phi0 = ParameterLib::findParameter<double>(
        //! \ogs_file_param{properties__property__VermaPruessModel__initial_porosity}
        config.getConfigParameter<std::string>("initial_porosity"), parameters,
        1, nullptr);

    auto const& phi_c = ParameterLib::findParameter<double>(
        //! \ogs_file_param{properties__property__VermaPruessModel__critical_porosity}
        config.getConfigParameter<std::string>("critical_porosity"), parameters,
        1, nullptr);

    auto const& n = ParameterLib::findParameter<double>(
        //! \ogs_file_param{properties__property__VermaPruessModel__exponent}
        config.getConfigParameter<std::string>("exponent"), parameters, 1,
        nullptr);

    return std::make_unique<VermaPruessModel>(k0, phi0, phi_c, n);
}
}  // namespace MaterialPropertyLib
