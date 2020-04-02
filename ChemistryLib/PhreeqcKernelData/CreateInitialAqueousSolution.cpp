/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateInitialAqueousSolution.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "InitialAqueousSolution.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
InitialAqueousSolution createInitialAqueousSolution(
    BaseLib::ConfigTree const& config,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map)
{
    std::map<std::string, cxxISolutionComp> components;
    //! \ogs_file_param{prj__chemical_system__solution__components}
    auto comp_config = config.getConfigSubtree("components");
    for (
        auto const& component_name :
        //! \ogs_file_param{prj__chemical_system__solution__components__component}
        comp_config.getConfigParameterList<std::string>("component"))
    {
        Component component(component_name);
        components.emplace(component_name, component);
    }

    for (auto const& component : components)
    {
        auto process_id_to_component_name_pair =
            std::find_if(process_id_to_component_name_map.begin(),
                         process_id_to_component_name_map.end(),
                         [&component](auto const& pair) {
                             return pair.second == component.first;
                         });

        if (process_id_to_component_name_pair ==
            process_id_to_component_name_map.end())
        {
            OGS_FATAL(
                "Component %s given in <solution>/<components> is not found in "
                "specified coupled processes (see "
                "<process>/<process_variables>/<concentration>).",
                component.first.c_str());
        }
    }

    if (components.size() + 1 != process_id_to_component_name_map.size())
    {
        OGS_FATAL(
            "The number of components given in <solution>/<components> is not "
            "in line with the number of transport processes - 1 which stands "
            "for the transport process of hydrogen.");
    }

    auto const means_of_adjusting_charge_in_str =
        //! \ogs_file_param{prj__chemical_system__solution__means_of_adjusting_charge}
        config.getConfigParameter<std::string>("means_of_adjusting_charge", "");

    if (means_of_adjusting_charge_in_str == "pH")
    {
        Component component("H(1)");
        component.Set_equation_name("charge");
        components.emplace("H(1)", component);
    }

    return InitialAqueousSolution(components);
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
