// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateSecondaryVariables.h"

#include "BaseLib/ConfigTree.h"
#include "SecondaryVariable.h"

namespace ProcessLib
{
void createSecondaryVariables(BaseLib::ConfigTree const& config,
                              SecondaryVariableCollection& secondary_variables)
{
    auto sec_vars_config =
        //! \ogs_file_param{prj__processes__process__secondary_variables}
        config.getConfigSubtreeOptional("secondary_variables");
    if (!sec_vars_config)
    {
        return;
    }

    for (
        auto sec_var_config :
        //! \ogs_file_param{prj__processes__process__secondary_variables__secondary_variable}
        sec_vars_config->getConfigSubtreeList("secondary_variable"))
    {
        auto const type =
            //! \ogs_file_attr{prj__processes__process__secondary_variables__secondary_variable__type}
            sec_var_config.getConfigAttributeOptional<std::string>("type");
        if (type)
        {
            WARN(
                "Secondary variable type specification is deprecated and is "
                "ignored. All secondary variable types are 'static'.");
        }

        auto const [internal_name, output_name] =
            [&sec_var_config]() -> std::pair<std::string, std::string>
        {
            auto const name =
                //! \ogs_file_attr{prj__processes__process__secondary_variables__secondary_variable__name}
                sec_var_config.getConfigAttributeOptional<std::string>("name");
            if (name)
            {
                return {*name, *name};
            }
            else
            {
                auto const internal_name =
                    //! \ogs_file_attr{prj__processes__process__secondary_variables__secondary_variable__internal_name}
                    sec_var_config.getConfigAttribute<std::string>(
                        "internal_name");
                auto const output_name =
                    //! \ogs_file_attr{prj__processes__process__secondary_variables__secondary_variable__output_name}
                    sec_var_config.getConfigAttribute<std::string>(
                        "output_name");
                return {internal_name, output_name};
            }
        }();

        secondary_variables.addNameMapping(internal_name, output_name);
    }
}

}  // namespace ProcessLib
