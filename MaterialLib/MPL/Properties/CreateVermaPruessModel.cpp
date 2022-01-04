/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
