/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/Utils.h"
#include "TemperatureDependentDiffusion.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createTemperatureDependentDiffusion(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "TemperatureDependentDiffusion");

    auto const& D0 = ParameterLib::findParameter<double>(
        //! \ogs_file_param{properties__property__TemperatureDependentDiffusion__reference_diffusion}
        config.getConfigParameter<std::string>("reference_diffusion"),
        parameters, 0, nullptr);

    auto const Ea =
        //! \ogs_file_param{properties__property__TemperatureDependentDiffusion__activation_energy}
        config.getConfigParameter<double>("activation_energy");

    auto const T0 =
        //! \ogs_file_param{properties__property__TemperatureDependentDiffusion__reference_temperature}
        config.getConfigParameter<double>("reference_temperature");

    return std::make_unique<TemperatureDependentDiffusion>(D0, Ea, T0);
}
}  // namespace MaterialPropertyLib
