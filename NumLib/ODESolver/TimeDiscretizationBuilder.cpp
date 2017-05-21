/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TimeDiscretizationBuilder.h"

#include "BaseLib/Error.h"

namespace NumLib
{
std::unique_ptr<TimeDiscretization> createTimeDiscretization(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__time_discretization__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "BackwardEuler")
    {
        //! \ogs_file_param_special{prj__time_loop__processes__process__time_discretization__BackwardEuler}
        return std::make_unique<BackwardEuler>();
    }
    else if (type == "ForwardEuler")
    {
        //! \ogs_file_param_special{prj__time_loop__processes__process__time_discretization__ForwardEuler}
        return std::make_unique<ForwardEuler>();
    }
    else if (type == "CrankNicolson")
    {
        //! \ogs_file_param{prj__time_loop__processes__process__time_discretization__CrankNicolson__theta}
        auto const theta = config.getConfigParameter<double>("theta");
        //! \ogs_file_param_special{prj__time_loop__processes__process__time_discretization__CrankNicolson}
        return std::make_unique<CrankNicolson>(theta);
    }
    else if (type == "BackwardDifferentiationFormula")
    {
        //! \ogs_file_param{prj__time_loop__processes__process__time_discretization__BackwardDifferentiationFormula__order}
        auto const order = config.getConfigParameter<unsigned>("order");
        //! \ogs_file_param_special{prj__time_loop__processes__process__time_discretization__BackwardDifferentiationFormula}
        return std::make_unique<BackwardDifferentiationFormula>(order);
    }
    else
    {
        OGS_FATAL("Unrecognized time discretization type `%s'", type.c_str());
    }
}
}
