/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSolutionComponent.h"
#include "AqueousSolution.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<Component> createSolutionComponents(
    BaseLib::ConfigTree const& config,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map)
{
    std::vector<Component> components;
    //! \ogs_file_param{prj__chemical_system__solution__components}
    auto comp_config = config.getConfigSubtree("components");
    for (
        auto const& component_name :
        //! \ogs_file_param{prj__chemical_system__solution__components__component}
        comp_config.getConfigParameterList<std::string>("component"))
    {
        components.emplace_back(component_name);
    }

    for (auto const& component : components)
    {
        auto process_id_to_component_name_map_element = std::find_if(
            process_id_to_component_name_map.begin(),
            process_id_to_component_name_map.end(),
            [&component](std::pair<int, std::string> const& map_element) {
                return map_element.second == component.name;
            });

        if (process_id_to_component_name_map_element ==
            process_id_to_component_name_map.end())
        {
            OGS_FATAL(
                "Component %s given in <solution>/<components> is not found in "
                "specified coupled processes (see "
                "<process>/<process_variables>/<concentration>).",
                component.name.c_str());
        }
    }
    if (components.size() + 1 != process_id_to_component_name_map.size())
    {
        OGS_FATAL(
            "The number of components given in <solution>/<components> is not "
            "in line with the number of transport processes - 1 which stands "
            "for the transport process of hydrogen.");
    }

    return components;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
