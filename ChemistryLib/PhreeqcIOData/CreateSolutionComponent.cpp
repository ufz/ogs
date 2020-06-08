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

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::vector<Component> createSolutionComponents(
    BaseLib::ConfigTree const& config, std::size_t const num_chemical_systems)
{
    std::vector<Component> components;
    //! \ogs_file_param{prj__chemical_system__solution__components}
    auto comp_config = config.getConfigSubtree("components");
    for (
        auto const& component_name :
        //! \ogs_file_param{prj__chemical_system__solution__components__component}
        comp_config.getConfigParameterList<std::string>("component"))
    {
        components.emplace_back(component_name, num_chemical_systems);
    }

    return components;
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
