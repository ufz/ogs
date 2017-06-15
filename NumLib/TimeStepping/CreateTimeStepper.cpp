/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateTimeStepper.cpp
 *  Created on May 2, 2017, 12:18 PM
 */

#include "CreateTimeStepper.h"

#include "BaseLib/ConfigTree.h"

#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"
#include "NumLib/TimeStepping/Algorithms/EvolutionaryPIDcontroller.h"

namespace NumLib
{
std::unique_ptr<TimeStepAlgorithm> createTimeStepper(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__time_loop__time_stepping__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    std::unique_ptr<NumLib::TimeStepAlgorithm> timestepper;

    if (type == "SingleStep")
    {
        //! \ogs_file_param_special{prj__time_loop__time_stepping__SingleStep}
        config.ignoreConfigParameter("type");
        timestepper.reset(new NumLib::FixedTimeStepping(0.0, 1.0, 1.0));
    }
    else if (type == "FixedTimeStepping")
    {
        timestepper = NumLib::FixedTimeStepping::newInstance(config);
    }
    else if (type == "EvolutionaryPIDcontroller")
    {
        timestepper = NumLib::createEvolutionaryPIDcontroller(config);
    }
    else
    {
        OGS_FATAL(
            "Unknown time stepping type: `%s'. "
            "The available types are \n\tSingleStep, \n\tFixedTimeStepping"
            "\n\tEvolutionaryPIDcontroller\n",
            type.data());
    }

    return timestepper;
}

}  // end of namespace NumLib