/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConstantParameter.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"

namespace ParameterLib
{
std::unique_ptr<ParameterBase> createConstantParameter(
    std::string const& name, BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "Constant");

    auto const value_data =
        //! \ogs_file_param{prj__parameters__parameter__Constant__value}
        config.getConfigParameterOptional<std::vector<double>>("value");

    auto const values_data =
        //! \ogs_file_param{prj__parameters__parameter__Constant__values}
        config.getConfigParameterOptional<std::vector<double>>("values");

    if ((value_data) && (values_data))
    {
        OGS_FATAL(
            "Both value and values tags were given.\n\
             Please give only one of them.");
    }

    std::vector<double> values;

    if (value_data)
    {
        values = *value_data;
    }
    else if (values_data)
    {
        values = *values_data;
    }

    if (values.empty())
    {
        OGS_FATAL("No value(s) was(were) provided.");
    }

    DBUG("Using following values for the constant parameter:");
    for (double const v : values)
    {
        (void)v;  // unused value if building w/o DBUG output.
        DBUG("\t{:g}", v);
    }

    return std::make_unique<ConstantParameter<double>>(name, std::move(values));
}

}  // namespace ParameterLib
