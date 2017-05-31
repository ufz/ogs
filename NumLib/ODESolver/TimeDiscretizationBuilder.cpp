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

    //! \ogs_file_param_special{prj__time_loop__processes__process__time_discretization__BackwardEuler}
    if (type == "BackwardEuler")
    {
        return std::make_unique<BackwardEuler>();
    }
    //! \ogs_file_param_special{prj__time_loop__processes__process__time_discretization__ForwardEuler}
    if (type == "ForwardEuler")
    {
        return std::make_unique<ForwardEuler>();
    }
    //! \ogs_file_param_special{prj__time_loop__processes__process__time_discretization__CrankNicolson}
    if (type == "CrankNicolson")
    {
        //! \ogs_file_param{prj__time_loop__processes__process__time_discretization__CrankNicolson__theta}
        auto const theta = config.getConfigParameter<double>("theta");
        return std::make_unique<CrankNicolson>(theta);
    }
    //! \ogs_file_param_special{prj__time_loop__processes__process__time_discretization__BackwardDifferentiationFormula}
    if (type == "BackwardDifferentiationFormula")
    {
        //! \ogs_file_param{prj__time_loop__processes__process__time_discretization__BackwardDifferentiationFormula__order}
        auto const order = config.getConfigParameter<unsigned>("order");
        return std::make_unique<BackwardDifferentiationFormula>(order);
    }

    OGS_FATAL("Unrecognized time discretization type `%s'", type.c_str());
}
}
