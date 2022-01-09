/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TimeDiscretizationBuilder.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "TimeDiscretization.h"

namespace NumLib
{
std::unique_ptr<TimeDiscretization> createTimeDiscretization(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__time_discretization__type}
    auto const type = config.getConfigParameter<std::string>("type");

    //! \ogs_file_param_special{prj__time_loop__processes__process__time_discretization__BackwardEuler}
    if (type == "BackwardEuler")
    {
        return std::make_unique<BackwardEuler>();
    }
    OGS_FATAL("Unrecognized time discretization type `{:s}'", type);
}
}  // namespace NumLib
