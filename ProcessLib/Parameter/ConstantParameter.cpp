/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConstantParameter.h"
#include <logog/include/logog.hpp>
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createConstantParameter(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{parameter__type}
    config.checkConfigParameter("type", "Constant");

    // Optional case for single-component variables where the value can be used.
    // If the 'value' tag is found, use it and return. Otherwise continue with
    // then required tag 'values'.
    {
        //! \ogs_file_param{parameter__Constant__value}
        auto const value = config.getConfigParameterOptional<double>("value");

        if (value)
        {
            DBUG("Using value %g for constant parameter.", *value);
            return std::unique_ptr<ParameterBase>(
                new ConstantParameter<double>(*value));
        }
    }

    // Value tag not available; continue with required values tag.
    std::vector<double> const values =
        //! \ogs_file_param{parameter__Constant__values}
        config.getConfigParameter<std::vector<double>>("values");

    if (values.empty())
        OGS_FATAL("No value available for constant parameter.");

    DBUG("Using following values for the constant parameter:");
    for (double const v : values)
    {
        (void)v;  // unused value if building w/o DBUG output.
        DBUG("\t%g", v);
    }

    return std::unique_ptr<ParameterBase>(
        new ConstantParameter<double>(values));
}

}  // ProcessLib
