/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on May 2, 2017, 12:18 PM
 */

#include "CreateTimeStepper.h"

#include <memory>
#include <string>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "NumLib/TimeStepping/Algorithms/CreateEvolutionaryPIDcontroller.h"
#include "NumLib/TimeStepping/Algorithms/CreateFixedTimeStepping.h"
#include "NumLib/TimeStepping/Algorithms/CreateIterationNumberBasedTimeStepping.h"
#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"

namespace NumLib
{
std::unique_ptr<TimeStepAlgorithm> createTimeStepper(
    BaseLib::ConfigTree const& config,
    std::vector<double> const& fixed_times_for_output)
{
    //! \ogs_file_param{prj__time_loop__processes__process__time_stepping__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "SingleStep")
    {
        //! \ogs_file_param_special{prj__time_loop__processes__process__time_stepping__SingleStep}
        config.ignoreConfigParameter("type");
        return std::make_unique<NumLib::FixedTimeStepping>(0.0, 1.0, 1.0);
    }
    if (type == "FixedTimeStepping")
    {
        return NumLib::createFixedTimeStepping(parseFixedTimeStepping(config),
                                               fixed_times_for_output);
    }
    if (type == "EvolutionaryPIDcontroller")
    {
        return NumLib::createEvolutionaryPIDcontroller(
            parseEvolutionaryPIDcontroller(config), fixed_times_for_output);
    }
    if (type == "IterationNumberBasedTimeStepping")
    {
        return NumLib::createIterationNumberBasedTimeStepping(
            parseIterationNumberBasedTimeStepping(config),
            fixed_times_for_output);
    }
    OGS_FATAL(
        "Unknown time stepping type: '{:s}'. The available types are: "
        "\n\tSingleStep,"
        "\n\tFixedTimeStepping,"
        "\n\tEvolutionaryPIDcontroller,",
        "\n\tIterationNumberBasedTimeStepping\n",
        type.data());
}

}  // end of namespace NumLib
