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

    // Optional case for single-component variables where the value can be used.
    // If the 'value' tag is found, use it and return. Otherwise continue with
    // then required tag 'values'.
    {
        auto const value =
            //! \ogs_file_param{prj__parameters__parameter__Constant__value}
            config.getConfigParameterOptional<std::vector<double>>("value");

        if (value)
        {
            if (value->size() != 1)
            {
                OGS_FATAL(
                    "Expected to read exactly one value, but {:d} were given.",
                    value->size());
            }
            DBUG("Using value {:g} for constant parameter.", (*value)[0]);
            return std::make_unique<ConstantParameter<double>>(name,
                                                               (*value)[0]);
        }
    }

    // Value tag not available; continue with required values tag.
    std::vector<double> const values =
        //! \ogs_file_param{prj__parameters__parameter__Constant__values}
        config.getConfigParameter<std::vector<double>>("values");

    if (values.empty())
    {
        OGS_FATAL("No value available for constant parameter.");
    }

    DBUG("Using following values for the constant parameter:");
    for (double const v : values)
    {
        (void)v;  // unused value if building w/o DBUG output.
        DBUG("\t{:g}", v);
    }

    return std::make_unique<ConstantParameter<double>>(name, values);
}

}  // namespace ParameterLib
