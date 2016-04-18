/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Process.h"

#include <algorithm>

#include <logog/include/logog.hpp>

#include "ProcessVariable.h"

namespace ProcessLib
{
ProcessVariable& findProcessVariable(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& pv_config, std::string const& tag)
{
	// Find process variable name in process config.
	std::string const name = pv_config.getConfParam<std::string>(tag);

        // Find corresponding variable by name.
	auto variable = std::find_if(variables.cbegin(), variables.cend(),
	                             [&name](ProcessVariable const& v)
	                             {
		                             return v.getName() == name;
	                             });

	if (variable == variables.end())
	{
		ERR(
		    "Could not find process variable '%s' in the provided variables "
		    "list for config tag <%s>.",
		    name.c_str(), tag.c_str());
		std::abort();
	}
	DBUG("Found process variable \'%s\' for config tag <%s>.",
		 variable->getName().c_str(), tag.c_str());

	// Const cast is needed because of variables argument constness.
	return const_cast<ProcessVariable&>(*variable);
}

std::vector<std::reference_wrapper<ProcessVariable>>
findProcessVariables(
		std::vector<ProcessVariable> const& variables,
		BaseLib::ConfigTree const& process_config,
		std::initializer_list<std::string> tag_names)
{
	std::vector<std::reference_wrapper<ProcessVariable>> vars;
	vars.reserve(tag_names.size());

	auto const pv_conf = process_config.getConfSubtree("process_variables");

	for (auto const& tag : tag_names) {
		vars.emplace_back(findProcessVariable(variables, pv_conf, tag));
	}

	return vars;
}

}  // namespace ProcessLib
