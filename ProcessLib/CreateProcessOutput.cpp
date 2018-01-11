/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateProcessOutput.h"

namespace ProcessLib
{
ProcessOutput createProcessOutput(BaseLib::ConfigTree const& output_config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__output__variables}
    auto const out_vars = output_config.getConfigSubtree("variables");

    std::set<std::string> output_variables;
    //! \ogs_file_param{prj__time_loop__processes__process__output__variables__variable}
    for (auto out_var :
         out_vars.getConfigParameterList<std::string>("variable"))
    {
        if (output_variables.find(out_var) != output_variables.cend())
        {
            OGS_FATAL("output variable `%s' specified more than once.",
                      out_var.c_str());
        }

        DBUG("adding output variable `%s'", out_var.c_str());
        output_variables.insert(out_var);
    }

    bool const output_residuals =
        //! \ogs_file_param{prj__time_loop__processes__process__output__output_extrapolation_residuals}
        output_config.getConfigParameter<bool>("output_extrapolation_residuals",
                                               false);

    return {output_variables, output_residuals};
}

}  // namespace ProcessLib
