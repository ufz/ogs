/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
        //! \ogs_file_param{properties__property__KozenyCarmanModel__initial_permeability}
        config.getConfigParameter<std::string>("initial_permeability"),
        parameters, 0, nullptr);

    auto const& phi0 = ParameterLib::findParameter<double>(
        //! \ogs_file_param{properties__property__KozenyCarmanModel__initial_porosity}
        config.getConfigParameter<std::string>("initial_porosity"), parameters,
        1, nullptr);

    return std::make_unique<KozenyCarmanModel>(k0, phi0);
}
}  // namespace MaterialPropertyLib
