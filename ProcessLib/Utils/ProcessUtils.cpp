/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessUtils.h"

#include "BaseLib/ConfigTree.h"
#include "ProcessLib/ProcessVariable.h"

namespace
{
ProcessLib::ProcessVariable& findVariableByName(
    std::vector<ProcessLib::ProcessVariable> const& variables,
    std::string const& name, std::string const& tag)
{
    // Find corresponding variable by name.
    auto variable = std::find_if(variables.cbegin(), variables.cend(),
                                 [&name](ProcessLib::ProcessVariable const& v)
                                 { return v.getName() == name; });

    if (variable == variables.end())
    {
        OGS_FATAL(
            "There is no entry of the defined process variable '{:s}' in the "
            "provided variables list (see tag <process_variables>). A "
            "definition is required for config tag <{:s}>.",
            name, tag);
    }
    DBUG("Found process variable '{:s}' for config tag <{:s}>.",
         variable->getName(), tag);

    // Const cast is needed because of variables argument constness.
    return const_cast<ProcessLib::ProcessVariable&>(*variable);
}
}  // namespace

namespace ProcessLib
{
ProcessVariable& findProcessVariable(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& pv_config, std::string const& tag)
{
    // Find process variable name in process config.
    //! \ogs_file_special
    std::string const name = pv_config.getConfigParameter<std::string>(tag);
    return findVariableByName(variables, name, tag);
}

std::vector<std::reference_wrapper<ProcessVariable>> findProcessVariables(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& pv_config,
    std::initializer_list<std::string>
        tags)
{
    std::vector<std::reference_wrapper<ProcessVariable>> vars;
    vars.reserve(variables.size());

    if (variables.size() > tags.size())
        DBUG("Found multiple process variables with a same tag.");

    for (auto& tag : tags)
    {
        auto vars_per_tag = findProcessVariables(variables, pv_config, tag);
        vars.insert(vars.end(), vars_per_tag.begin(), vars_per_tag.end());
    }

    return vars;
}

std::vector<std::reference_wrapper<ProcessVariable>> findProcessVariables(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& pv_config,
    std::string const& tag)
{
    std::vector<std::reference_wrapper<ProcessVariable>> vars;

    //! \ogs_file_special
    auto var_names = pv_config.getConfigParameterList<std::string>(tag);

    if (var_names.empty())
    {
        OGS_FATAL("No entity is found with config tag <{:s}>.", tag);
    }

    std::vector<std::string> cached_var_names;

    for (std::string const& var_name : var_names)
    {
        vars.emplace_back(findVariableByName(variables, var_name, tag));
        cached_var_names.push_back(var_name);
    }

    // Eliminate duplicates in the set of variable names
    BaseLib::makeVectorUnique(cached_var_names);

    if (cached_var_names.size() != var_names.size())
    {
        OGS_FATAL("Found duplicates with config tag <{:s}>.", tag);
    }

    return vars;
}
}  // namespace ProcessLib
