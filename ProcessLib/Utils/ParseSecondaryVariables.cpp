/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ParseSecondaryVariables.h"
#include "BaseLib/ConfigTree.h"
#include "NumLib/NamedFunctionCaller.h"
#include "ProcessLib/SecondaryVariable.h"

namespace ProcessLib
{
void parseSecondaryVariables(
    BaseLib::ConfigTree const& config,
    SecondaryVariableCollection& secondary_variables,
    NumLib::NamedFunctionCaller& named_function_caller)
{
    auto sec_vars_config =
        config.getConfigSubtreeOptional("secondary_variables");
    if (!sec_vars_config)
        return;

    for (auto sec_var_config :
         sec_vars_config->getConfigSubtreeList("secondary_variable")) {
        auto const type =
            sec_var_config.getConfigAttribute<std::string>("type");

        auto const internal_name =
            sec_var_config.getConfigAttribute<std::string>("internal");
        auto const external_name =
            sec_var_config.getConfigAttribute<std::string>("external");

        secondary_variables.addNameMapping(internal_name, external_name);

        if (type == "static") {
            // alright
        } else if (type == "dynamic") {
            auto const& sink_fct = internal_name;

            for (auto const plug :
                 sec_var_config.getConfigParameterList("plug")) {
                auto const sink_arg =
                    plug.getConfigAttribute<std::string>("sink_arg");
                auto const source =
                    plug.getConfigAttribute<std::string>("source");

                named_function_caller.plug(sink_fct, sink_arg, source);
            }
        }
    }
}

}  // namespace ProcessLib
