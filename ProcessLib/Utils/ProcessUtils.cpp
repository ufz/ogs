/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessUtils.h"
#include <iterator>
#include "ProcessLib/ProcessVariable.h"

namespace ProcessLib
{
ProcessVariable& findProcessVariable(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& pv_config, std::string const& tag)
{
    // Find process variable name in process config.
    //! \ogs_file_special
    std::string const name = pv_config.getConfigParameter<std::string>(tag);

    // Find corresponding variable by name.
    auto variable = std::find_if(
        variables.cbegin(), variables.cend(),
        [&name](ProcessVariable const& v) { return v.getName() == name; });

    if (variable == variables.end())
    {
        OGS_FATAL(
            "Could not find process variable '%s' in the provided variables "
            "list for config tag <%s>.",
            name.c_str(), tag.c_str());
    }
    DBUG("Found process variable '%s' for config tag <%s>.",
         variable->getName().c_str(), tag.c_str());

    // Const cast is needed because of variables argument constness.
    return const_cast<ProcessVariable&>(*variable);
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

    auto var_names = pv_config.getConfigParameterList<std::string>(tag);

    if (var_names.size() == 0)
        OGS_FATAL("Config tag <%s> is not found.", tag.c_str());

    // collect variable names to check if there are duplicates
    std::set<std::string> cached_var_names;

    for (std::string const& var_name : var_names)
    {
        auto variable = std::find_if(variables.cbegin(), variables.cend(),
                                     [&var_name](ProcessVariable const& v) {
                                         return v.getName() == var_name;
                                     });

        if (variable == variables.end())
        {
            OGS_FATAL(
                "Could not find process variable '%s' in the provided "
                "variables "
                "list for config tag <%s>.",
                var_name.c_str(), tag.c_str());
        }
        DBUG("Found process variable \'%s\' for config tag <%s>.",
             variable->getName().c_str(), tag.c_str());

        vars.emplace_back(const_cast<ProcessVariable&>(*variable));

        cached_var_names.insert(var_name);
    }

    if (cached_var_names.size() != var_names.size())
        OGS_FATAL("Found duplicates with config tag <%s>.", tag.c_str());

    return vars;
}
}  // ProcessLib
