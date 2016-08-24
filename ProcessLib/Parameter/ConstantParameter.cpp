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

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createConstantParameter(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{parameter__type}
    config.checkConfigParameter("type", "Constant");
    //! \ogs_file_param{parameter__Constant__value}
    auto const value = config.getConfigParameter<double>("value");
    DBUG("Using value %g", value);

    return std::unique_ptr<ParameterBase>(new ConstantParameter<double>(value));
}

}  // ProcessLib
