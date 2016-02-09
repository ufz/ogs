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
    BaseLib::ConfigTree const& process_config, std::string const& tag,
    std::vector<ProcessVariable> const& variables)
{
	// Find process variable name in process config.
	std::string const name = process_config.getConfParam<std::string>(tag);

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
	DBUG("Found process variable \'%s\'.", variable->getName().c_str());

	// Const cast is needed because of variables argument constness.
	return const_cast<ProcessVariable&>(*variable);
}

}  // namespace ProcessLib
