/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 10, 2019
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "Constant.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Constant> createConstant(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Constant");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create Constant property {:s}.", property_name);
    auto const value_data =
        //! \ogs_file_param{properties__property__Constant__value}
        config.getConfigParameterOptional<std::vector<double>>("value");

    auto const values_data =
        //! \ogs_file_param{properties__property__Constant__values}
        config.getConfigParameterOptional<std::vector<double>>("values");

    if ((value_data) && (values_data))
    {
        OGS_FATAL(
            "Both the value and the values tags were given.\n\
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

    return std::make_unique<Constant>(std::move(property_name),
                                      fromVector(values));
}
}  // namespace MaterialPropertyLib
