// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateSolutionComponent.h"

#include "AqueousSolution.h"
#include "BaseLib/ConfigTree.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<Component> createSolutionComponents(
    BaseLib::ConfigTree const& config)
{
    std::vector<Component> components;
    //! \ogs_file_param{prj__chemical_system__solution__components}
    auto components_config = config.getConfigSubtree("components");

    for (
        auto const& component_config :
        //! \ogs_file_param{prj__chemical_system__solution__components__component}
        components_config.getConfigSubtreeList("component"))
    {
        auto const component_name = component_config.getValue<std::string>();
        auto const chemical_formula =
            //! \ogs_file_attr{prj__chemical_system__solution__components__component__chemical_formula}
            component_config.getConfigAttribute<std::string>("chemical_formula",
                                                             "");
        components.emplace_back(component_name, chemical_formula);
    }

    return components;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
