/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   CreateComponentProperties.cpp
 */

#include "CreateComponentProperties.h"
#include "BaseLib/ConfigTree.h"
#include "ComponentProperties.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace MaterialLib
{
namespace Component
{
std::vector<ComponentProperties> createComponentProperties(
    BaseLib::ConfigTree const& configs,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    DBUG("Create PerComponentProperties.");
    std::vector<ComponentProperties> per_component_properties;

    auto const& comp_configs =
        //! \ogs_file_param{material__components}
        configs.getConfigSubtree("components");

    for (auto const& comp_config :
         //! \ogs_file_param{material__components__component}
         comp_configs.getConfigSubtreeList("component"))
    {
        //! \ogs_file_param{material__components__component__name}
        auto const name = comp_config.getConfigAttribute<std::string>("name");

        auto const& Dm = ProcessLib::findParameter<double>(
            comp_config,
            //! \ogs_file_param_special{material__component__molecular_diffusion}
            "molecular_diffusion", parameters, 1);

        auto const& beta_l = ProcessLib::findParameter<double>(
            comp_config,
            //! \ogs_file_param_special{material__component__solute_dispersivity_longitudinal}
            "solute_dispersivity_longitudinal", parameters, 1);

        auto const& beta_t = ProcessLib::findParameter<double>(
            comp_config,
            //! \ogs_file_param_special{material__component__solute_dispersivity_transverse}
            "solute_dispersivity_transverse", parameters, 1);

        per_component_properties.push_back(
            ComponentProperties{name, Dm, beta_l, beta_t});
    }

    return per_component_properties;
}
}  // namespace Component
}  // namespace MaterialLib
