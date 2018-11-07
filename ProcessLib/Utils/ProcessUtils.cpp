/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessUtils.h"
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
    DBUG("Found process variable \'%s\' for config tag <%s>.",
         variable->getName().c_str(), tag.c_str());

    // Const cast is needed because of variables argument constness.
    return const_cast<ProcessVariable&>(*variable);
}

std::vector<std::reference_wrapper<ProcessVariable>> findProcessVariables(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& pv_config,
    std::initializer_list<std::string>
        tag_names)
{
    std::vector<std::reference_wrapper<ProcessVariable>> vars;
    vars.reserve(variables.size());

    if (variables.size() > tag_names.size())
    {
        DBUG("Found %d process variables with only %d config tags.",
             variables.size(), tag_names.size());

        for (auto const& tag : tag_names)
        {
            auto vars_per_tag = findProcessVariables(variables, pv_config, tag);

            vars.insert(vars.end(), vars_per_tag.begin(), vars_per_tag.end());
        }
    }
    else
    {
        for (auto const& tag : tag_names)
        {
            vars.emplace_back(findProcessVariable(variables, pv_config, tag));
        }
    }

    return vars;
}

std::vector<std::reference_wrapper<ProcessVariable>> findProcessVariables(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& pv_config,
    std::string const& tag_name)
{
    std::vector<std::reference_wrapper<ProcessVariable>> vars;

    auto range = pv_config.getConfigParameterList<std::string>(tag_name);

    for (std::string const& name : range)
    {
        auto variable = std::find_if(
            variables.cbegin(), variables.cend(),
            [&name](ProcessVariable const& v) { return v.getName() == name; });

        if (variable == variables.end())
        {
            OGS_FATAL(
                "Could not find process variable '%s' in the provided variables "
                "list for config tag <%s>.",
                name.c_str(), tag_name.c_str());
        }
        DBUG("Found process variable \'%s\' for config tag <%s>.",
             variable->getName().c_str(), tag_name.c_str());

        vars.emplace_back(const_cast<ProcessVariable&>(*variable));
    }

    return vars;

}
}  // ProcessLib
